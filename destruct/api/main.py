import argparse

import destruct.api.run
import destruct.api.create_ref_data


def main():
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers()
    
    destruct.api.run.add_arguments(subparsers.add_parser('run'))
    destruct.api.create_ref_data.add_arguments(subparsers.add_parser('create_ref_data'))

    args = vars(argparser.parse_args())
    func = args.pop('func')
    func(**args)
