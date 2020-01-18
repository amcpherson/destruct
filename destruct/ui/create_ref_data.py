import argparse
import destruct.create_ref_data


def create_ref_data(**args):
    config = {}
    if args['config'] is not None:
        exec(open(args['config']).read(), {}, config)

    destruct.create_ref_data.create_ref_data(config, args['ref_data_dir'])


def add_arguments(argparser):
    argparser.add_argument('ref_data_dir',
                           help='Reference dataset directory')

    argparser.add_argument('-c', '--config',
                           help='Configuration filename')

    argparser.set_defaults(func=create_ref_data)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()

    add_arguments(argparser)

    args = vars(argparser.parse_args())
    func = args.pop('func')
    func(**args)
