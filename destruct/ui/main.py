import argparse

import destruct.ui.run
import destruct.ui.create_ref_data
import destruct.ui.extract_somatic


def main():
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers()
    
    destruct.ui.run.add_arguments(subparsers.add_parser('run'))
    destruct.ui.create_ref_data.add_arguments(subparsers.add_parser('create_ref_data'))
    destruct.ui.extract_somatic.add_arguments(subparsers.add_parser('extract_somatic'))

    args = vars(argparser.parse_args())
    func = args.pop('func')
    func(**args)
