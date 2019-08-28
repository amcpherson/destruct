import pandas as pd

import destruct.balanced


def classify_rearrangement_type(entry):
    break1 = entry['position_1']
    break2 = entry['position_2']
    size = abs(int(break1) - int(break2))
    orientation_type = entry['type']
    if entry['balanced']:
        rearrangement_type = 'balanced'
    elif size <= 1000000 and orientation_type == 'deletion':
        rearrangement_type = 'deletion'
    elif size <= 10000 and orientation_type == 'inversion':
        rearrangement_type = 'foldback'
    elif size <= 1000000 and orientation_type == 'inversion':
        rearrangement_type = 'inversion'
    elif size <= 1000000 and orientation_type == 'duplication':
        rearrangement_type = 'duplication'
    else:
        rearrangement_type = 'unbalanced'
    return rearrangement_type


def annotate_rearrangement_type(df):
    rearrangement_types = list()
    for idx, row in df.iterrows():
        rearrangement_types.append(classify_rearrangement_type(row))
    df['rearrangement_type'] = rearrangement_types


def filter_annotate_breakpoints(
        input_breakpoint_filename,
        input_breakpoint_library_filename,
        control_ids,
        output_breakpoint_filename,
        output_breakpoint_library_filename,
        patient_library_filename=None):
    """ Filter and annotate breakpoints.

    Args:
        input_breakpoint_filename (str): filename of breakpoint table
        input_breakpoint_library_filename (str): filename of breakpoint library table
        control_ids (list): control id or ids
        output_breakpoint_filename (str): output filename of breakpoint table
        output_breakpoint_library_filename (str): output filename of breakpoint library table

    KwArgs:
        patient_libraries (str): dataframe of library ids for each patient, columns patient_id, library

    If patient_libraries is not specified, assumed one patient for all libraries.

    """

    brk = pd.read_csv(input_breakpoint_filename, sep='\t',
                      converters={'chromosome_1': str, 'chromosome_2': str})

    brklib = pd.read_csv(input_breakpoint_library_filename, sep='\t')

    if patient_library_filename is None:
        patient_libraries = {'null': list(brklib['library'].unique())}

    else:
        patient_library_data = pd.read_csv(patient_library_filename, sep='\t')
        patient_libraries = {}
        for patient_id, rows in patient_library_data.groupby('patient_id'):
            patient_libraries[patient_id] = [row['library'] for row in rows]

    # Add is_normal column
    brklib['is_normal'] = brklib['library'].isin(control_ids)

    # Add patient_id column
    brklib['patient_id'] = ''
    for patient_id, library_ids in patient_libraries.items():
        for library_id in library_ids:
            brklib.loc[brklib['library'] == library_id, 'patient_id'] = patient_id

    num_patients = (
        brklib[['prediction_id', 'patient_id']]
        .drop_duplicates()
        .groupby('prediction_id')
        .size()
        .rename('num_patients')
    )

    # Mark as germline any prediction with nonzero normal reads
    is_germline = brklib[brklib['num_reads'] > 0].groupby('prediction_id')['is_normal'].any()

    # DGV predictions also germline
    is_dgv = brk.set_index('prediction_id')['dgv_ids'].notnull()

    brk.set_index('prediction_id', inplace=True)
    brk['is_germline'] = is_germline
    brk['is_dgv'] = is_dgv
    brk['num_patients'] = num_patients
    brk.reset_index(inplace=True)

    # Mark as a filtered any breakpoint that is common or germline
    brk['is_filtered'] = brk['is_germline'] | brk['is_dgv'] | (brk['num_patients'] > 1)

    # Get a list of breakends and their filtered status
    def get_brkend(brk, side, data_cols):
        cols = ['chromosome', 'strand', 'position']
        side_cols = [a + '_' + side for a in cols]
        brkend = brk[['prediction_id'] + side_cols + data_cols]
        brkend = brkend.rename(columns=dict(zip(side_cols, cols)))
        brkend['side'] = side
        return brkend

    brkend = pd.concat([get_brkend(brk, '1', ['is_filtered']),
                        get_brkend(brk, '2', ['is_filtered'])], ignore_index=True)

    # Calculate the minimum distance to the nearest filtered breakend for all breakpoints
    def calculate_dist_filtered(data):

        data = data.sort_values('position')

        data['nearest_left'] = None
        data['nearest_right'] = None

        data.loc[data['is_filtered'], 'nearest_left'] = data.loc[data['is_filtered'], 'position']
        data.loc[data['is_filtered'], 'nearest_right'] = data.loc[data['is_filtered'], 'position']

        data['nearest_left'] = data['nearest_left'].fillna(method='ffill')
        data['nearest_right'] = data['nearest_right'].fillna(method='bfill')

        data['nearest_left'] = data['position'] - data['nearest_left']
        data['nearest_right'] = data['nearest_right'] - data['position']

        data['dist_filtered'] = data[['nearest_left', 'nearest_right']].min(axis=1)

        return data['dist_filtered']

    dist_filtered = brkend.set_index(['prediction_id', 'side'])\
                          .groupby(['chromosome', 'strand'])\
                          .apply(calculate_dist_filtered)\
                          .reset_index(level=[0, 1], drop=True)\
                          .unstack(level=1)\
                          .min(axis=1)

    brk.set_index('prediction_id', inplace=True)
    brk['dist_filtered'] = dist_filtered
    brk.reset_index(inplace=True)

    chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                   '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
                   'X', 'Y']

    # Filter breakpoints
    brk = brk.loc[(~brk['is_filtered']) &
                  (brk['dist_filtered'] > 50) &
                  (brk['num_reads'] >= 2) &
                  (brk['num_split'] >= 2) &
                  (brk['log_likelihood'] > -20.) &
                  (brk['template_length_min'] > 120) &
                  (brk['chromosome_1'].isin(chromosomes)) &
                  (brk['chromosome_2'].isin(chromosomes))]

    brklib = brklib.merge(brk[['prediction_id']].drop_duplicates())

    # Balanced rearrangement annotation
    balanced_prediction_ids = []

    for patient_id, library_ids in patient_libraries.items():
        patient_prediction_ids = brklib.loc[brklib['library'].isin(library_ids), ['prediction_id']].drop_duplicates()
        patient_brks = brk.merge(patient_prediction_ids)
        balanced_rearrangements = destruct.balanced.detect_balanced_rearrangements(patient_brks)
        for rearrangement in balanced_rearrangements:
            balanced_prediction_ids.extend(rearrangement.prediction_ids)

    brk['balanced'] = False
    brk.loc[brk['prediction_id'].isin(balanced_prediction_ids), 'balanced'] = True

    # Annotate rearrangement type
    annotate_rearrangement_type(brk)

    # Output filtered breakpoint tables
    brk.to_csv(output_breakpoint_filename, sep='\t', header=True, index=False)
    brklib.to_csv(output_breakpoint_library_filename, sep='\t', header=True, index=False)


