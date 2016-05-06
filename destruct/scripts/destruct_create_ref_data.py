import argparse
import destruct.create_ref_data


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()

    argparser.add_argument('ref_data_dir',
                           help='Reference dataset directory')

    argparser.add_argument('-c', '--config',
                           help='Configuration filename')

    args = argparser.parse_args()

    args = vars(argparser.parse_args())

    config = {}
    if args['config'] is not None:
        execfile(args['config'], {}, config)

    destruct.create_ref_data.create_ref_data(config, args['ref_data_dir'])


