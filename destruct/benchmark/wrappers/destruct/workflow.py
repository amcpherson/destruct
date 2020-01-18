import pypeliner
import pypeliner.managed as mgd

import destruct.workflow
import destruct.create_ref_data
import destruct.benchmark.wrappers.destruct.tasks


def create_destruct_wrapper_workflow(bam_filenames, output_filename, raw_data_dir, control_id=None, config=None, ref_data_dir=None):
    workflow = pypeliner.workflow.Workflow(default_ctx={'mem': 4})

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=list(bam_filenames.keys()),
    )

    workflow.subworkflow(
        name='run_destruct',
        func=destruct.workflow.create_destruct_workflow,
        args=(
            mgd.InputFile('bam', 'sample_id', fnames=bam_filenames),
            mgd.TempOutputFile('breakpoint_table'),
            mgd.TempOutputFile('breakpoint_library_table'),
            mgd.TempOutputFile('breakpoint_read_table'),
            config,
            ref_data_dir,
        ),
        kwargs={
            'raw_data_dir': raw_data_dir,
        },
    )

    workflow.transform(
        name='post_process',
        func=destruct.benchmark.wrappers.destruct.tasks.destruct_postprocess,
        args=(
            mgd.TempInputFile('breakpoint_table'),
            mgd.TempInputFile('breakpoint_library_table'),
            mgd.OutputFile(output_filename),
        ),
        kwargs={
            'control_id': control_id,
        }
    )

    return workflow


def setup_destruct(test_config, **kwargs):
    config = kwargs.get('config', {})
    
    if test_config.get('chromosomes', None) is not None:
        config = {
            'chromosomes': test_config['chromosomes'],
            'ensembl_assemblies': ['chromosome.{}'.format(c) for c in test_config['chromosomes']],
        }

    destruct.create_ref_data.create_ref_data(config, kwargs['ref_data_dir'])

