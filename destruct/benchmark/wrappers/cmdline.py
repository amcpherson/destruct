import argparse


def interface(Wrapper):

    parser = argparse.ArgumentParser()

    parser.add_argument('install_directory', help='Installation directory')

    subparsers = parser.add_subparsers(dest='command')

    install_parser = subparsers.add_parser('install')
    
    install_parser.add_argument('--chromosomes', nargs='*', type=str, default=None,
                                help='Reference chromosomes')

    run_parser = subparsers.add_parser('run')

    run_parser.add_argument('output_filename',
                            help='Results output filename')

    run_parser.add_argument('temp_directory',
                            help='Temporary directory')

    run_parser.add_argument('--bam_files', nargs='+', required=True,
                            help='Input bam filenames')

    run_parser.add_argument('--lib_ids', nargs='+', required=True,
                            help='Input ids for respective bam filenames')

    args = parser.parse_args()

    wrapper = Wrapper(args.install_directory)

    if args.command == 'install':

        wrapper.install(chromosomes=args.chromosomes)

    else:

        if len(args['bam_files']) != len(args['lib_ids']):
            raise Exception('--lib_ids must correspond one to one with --bam_files')

        bam_filenames = dict(zip(args['lib_ids'], args['bam_files']))

        wrapper.run(bams, args.output_filename, args.temp_directory)





