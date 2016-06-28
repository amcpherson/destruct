import argparse

import destruct.ui.run
import destruct.ui.create_ref_data


def main():
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers()
    
    destruct.ui.run.add_arguments(subparsers.add_parser('run'))
    destruct.ui.create_ref_data.add_arguments(subparsers.add_parser('create_ref_data'))

    args = vars(argparser.parse_args())
    func = args.pop('func')
    func(**args)
