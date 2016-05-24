import pandas as pd


def destruct_postprocess(breakpoint_table_filename, breakpoint_library_table_filename, output_filename, control_id=None):

    breakpoint_table = pd.read_csv(breakpoint_table_filename, sep='\t')
    breakpoint_library_table = pd.read_csv(breakpoint_library_table_filename, sep='\t')

    breakpoint_counts = (
        breakpoint_library_table
        .set_index(['prediction_id', 'library'])[['num_reads']]
        .unstack(fill_value=0)
    )

    breakpoint_counts.columns = [a[1] + '_count' for a in breakpoint_counts.columns]

    breakpoint_table = breakpoint_table.merge(breakpoint_counts,
                                              left_on='prediction_id',
                                              right_index=True)
    
    # Filter based on evidence in control dataset
    if control_id is not None:
         breakpoint_table = breakpoint_table[breakpoint_table['{0}_count'.format(control_id)] == 0]          

    breakpoint_table.to_csv(output_filename, sep='\t', index=False)

