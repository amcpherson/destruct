import argparse


def interface(Wrapper):

    parser = argparse.ArgumentParser()

    parser.add_argument('install_directory', help='Installation directory')

    subparsers = parser.add_subparsers(dest='command')

    install_parser = subparsers.add_parser('install')
    install_parser.add_argument('--chromosomes', nargs='*', type=str, default=None, help='Reference chromosomes')

    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('tumour_bam', help='Tumour bam filename')
    run_parser.add_argument('normal_bam', help='Normal bam filename')
    run_parser.add_argument('output_filename', help='Results output filename')
    run_parser.add_argument('temp_directory', help='Temporary directory')

    args = parser.parse_args()

    wrapper = Wrapper(args.install_directory)

    if args.command == 'install':

        wrapper.install(chromosomes=args.chromosomes)

    else:

        wrapper.run(args.tumour_bam, args.normal_bam, args.output_filename, args.temp_directory)





