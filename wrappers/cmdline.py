import argparse


def interface(Wrapper):

    parser = argparse.ArgumentParser()

    parser.add_argument('install_directory', help='Installation directory')

    subparsers = parser.add_subparsers(dest='command')

    install_parser = subparsers.add_parser('install')
    install_parser.add_argument('--chromosomes', nargs='*', type=str, default=None, help='Reference chromosomes')

    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('temp_directory', help='Temporary directory')
    run_parser.add_argument('output_filename', help='Results output filename')
    run_parser.add_argument('--bam_filenames', required=True, nargs='+', type=str, default=None, help='Bam filenames')
    run_parser.add_argument('--sample_ids', required=True, nargs='+', type=str, default=None, help='Sample identifiers')

    args = parser.parse_args()

    wrapper = Wrapper(args.install_directory)

    if args.command == 'install':

        wrapper.install(chromosomes=args.chromosomes)

    else:

        if len(args.sample_ids) != len(args.bam_filenames):
            raise Exception('require one sample id per bam')

        wrapper.run(args.temp_directory, dict(zip(args.sample_ids, args.bam_filenames)), args.output_filename)





