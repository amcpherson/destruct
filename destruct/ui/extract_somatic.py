import argparse
import destruct.somatic


def filter_annotate_breakpoints(**args):
    destruct.somatic.filter_annotate_breakpoints(
        args['input_breakpoint_filename'],
        args['input_breakpoint_library_filename'],
        args['control_ids'],
        args['output_breakpoint_filename'],
        args['output_breakpoint_library_filename'],
        patient_library_filename=args['patient_libraries_filename'],
    )


def add_arguments(argparser):
    argparser.add_argument('input_breakpoint_filename',
                           help='Input breakpoints filename')

    argparser.add_argument('input_breakpoint_library_filename',
                           help='Input breakpoint library filename')

    argparser.add_argument('output_breakpoint_filename',
                           help='Output breakpoints filename')

    argparser.add_argument('output_breakpoint_library_filename',
                           help='Output breakpoint library filename')

    argparser.add_argument('--control_ids', nargs='+',
                           help='Control/normal library ids')

    argparser.add_argument('--patient_libraries_filename', required=False,
                           help='Mapping from patient_id to library as tsv')

    argparser.set_defaults(func=filter_annotate_breakpoints)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()

    add_arguments(argparser)

    args = vars(argparser.parse_args())
    func = args.pop('func')
    func(**args)
