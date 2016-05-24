import numpy as np

import destruct.balanced


def annotate_results(brks):
    # Classify balanced rearrangements as inversions, reciprocal translocations, complex
    balanced_rearrangements = destruct.balanced.detect_balanced_rearrangements(brks)

    brks.set_index('prediction_id', inplace=True)
    brks['rearrangement_type'] = ''

    for rearrangement in balanced_rearrangements:
        pred_types = brks.loc[rearrangement.prediction_ids, 'type'].values
        
        if len(rearrangement.prediction_ids) == 2 and np.all(pred_types == 'inversion'):
            brks.loc[rearrangement.prediction_ids, 'rearrangement_type'] == 'inversion'
            
        elif len(rearrangement.prediction_ids) == 2 and np.all(pred_types == 'translocation'):
            brks.loc[rearrangement.prediction_ids, 'rearrangement_type'] == 'reciprocal'

        else:
            brks.loc[rearrangement.prediction_ids, 'rearrangement_type'] == 'complex'

    brks.reset_index(inplace=True)

    # Classify non inversions with small distance as foldbacks
    brks['dist'] = np.absolute(brks['position_1'] - brks['position_2'])
    is_foldback = (brks['dist'] < 1000) & (brks['type'] == 'inversion') & (brks['rearrangement_type'] != 'inversion')
    brks.loc[is_foldback, 'rearrangement_type'] = 'foldback'

    # All non-complex deletions as deletions
    is_deletion = (brks['type'] == 'deletion') & (brks['rearrangement_type'] != 'complex')
    brks.loc[is_deletion, 'rearrangement_type'] = 'deletion'

    # All non-complex duplications as duplications
    is_duplication = (brks['type'] == 'duplication') & (brks['rearrangement_type'] != 'complex')
    brks.loc[is_duplication, 'rearrangement_type'] = 'duplication'

    # All non-complex non reciprocal translocations as translocations
    is_translocation = (brks['type'] == 'translocation') & (brks['rearrangement_type'] != 'complex') & (brks['rearrangement_type'] != 'reciprocal')
    brks.loc[is_translocation, 'rearrangement_type'] = 'translocation'

    # All non-complex non balanced inversions as unbalanced
    is_unbalanced = (brks['type'] == 'inversion') & (brks['rearrangement_type'] != 'complex') & (brks['rearrangement_type'] != 'inversion') & (brks['rearrangement_type'] != 'foldback')
    brks.loc[is_unbalanced, 'rearrangement_type'] = 'unbalanced'

    return brks

