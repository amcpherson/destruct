import pandas as pd

def read_select_write(df_iter, select, out_filename):
    """ Write subset of streamed data.

    Args:
        df_iter (iter of pandas.DataFrame): streamed data
        select (pandas.DataFrame): subset for filtering
        out_filename (str): output tsv filename

    """

    with open(out_filename, 'w') as out_file:
        for chunk in df_iter:
            chunk = chunk.merge(select, how='inner')
            chunk.to_csv(out_file, sep='\t', header=False, index=False)


def read_filter_write(df_iter, f_filter, out_filename):
    """ Write filtered set of streamed data.

    Args:
        df_iter (iter of pandas.DataFrame): streamed data
        f_filter (callable): filtering function
        out_filename (str): output tsv filename

    """

    with open(out_filename, 'w') as out_file:
        for chunk in df_iter:
            chunk = f_filter(chunk)
            chunk.to_csv(out_file, sep='\t', header=False, index=False)


def group_aware_iter(df_iter, group_cols):
    """ Stream data on group boundaries.

    Args:
        df_iter (iter of pandas.DataFrame): streamed data
        group_cols (list of str): columns defining group

    Yields:
        pandas.DataFrame: group aware data

    """

    prev_data = None

    for df in df_iter:

        if len(df.index) == 0:
            return

        last_group = df.loc[df.index[-1], group_cols].values

        # Remove last group
        next_data = df.loc[(df[group_cols] != last_group).all(axis=1)]

        # Add previous last group to beginning of data
        if prev_data is not None:
            next_data = pd.concat([prev_data, next_data], ignore_index=True)

        yield next_data

        # Save last group for next yield
        prev_data = df.loc[(df[group_cols] == last_group).all(axis=1)]

    yield prev_data


