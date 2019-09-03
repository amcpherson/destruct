import argparse
import yaml

import destruct.benchmark.align.bwa.workflow
import destruct.benchmark.destruct_test
import destruct.benchmark.create_breakpoint_simulation


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()

    argparser.add_argument('ref_data_dir',
                           help='Reference genomes directory')

    argparser.add_argument('tool_defs',
                           help='Tool Definition Filename')

    argparser.add_argument('--chromosomes', nargs='*', type=str, default=None,
                           help='Reference chromosomes')

    args = vars(argparser.parse_args())

    yaml_text = open(args['tool_defs']).read().format(ref_data_dir=args['ref_data_dir'])
    tool_defs = yaml.load(yaml_text)

    ctx = {'mem':4}

    test_config = {}
    if args.get('chromosomes', None) is not None:
        test_config['chromosomes'] = args['chromosomes']

    for tool_name, tool_info in tool_defs.items():
        destruct.benchmark.destruct_test.run_setup_function(tool_info, test_config)

