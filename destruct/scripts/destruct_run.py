import argparse

import pypeliner

import destruct
import destruct.workflow


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()

    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('ref_data_dir',
                           help='Reference dataset directory')

    argparser.add_argument('breakpoint_table',
                           help='Output table of breakpoint information in TSV format')

    argparser.add_argument('breakpoint_library_table',
                           help='Output table of library specific breakpoint information in TSV format')

    argparser.add_argument('--breakpoint_read_table', required=False, default=None,
                           help='Output table of breakpoint read information in TSV format')

    argparser.add_argument('--bam_files', nargs='+',
                           help='Input bam filenames')

    argparser.add_argument('--lib_ids', nargs='+',
                           help='Input ids for respective bam filenames')

    argparser.add_argument('--config', required=False,
                           help='Configuration filename')

    args = vars(argparser.parse_args())

    if len(args['bam_files']) != len(args['lib_ids']):
        raise Exception('--lib_ids must correspond one to one with --bam_files')

    bam_filenames = dict(zip(args['lib_ids'], args['bam_files']))

    config = {}
    if args['config'] is not None:
        execfile(args['config'], {}, config)

    pyp = pypeliner.app.Pypeline(modules=[destruct], config=args)

    workflow = destruct.workflow.create_destruct_workflow(
        bam_filenames,
        args['breakpoint_table'],
        args['breakpoint_library_table'],
        config,
        args['ref_data_dir'],
        breakpoint_read_table=args['breakpoint_read_table'],
    )

    pyp.run(workflow)


