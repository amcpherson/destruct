import os
import argparse
import yaml

import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

import destruct.benchmark.align.bwa.workflow
import destruct.benchmark.destruct_test
import destruct.benchmark.create_breakpoint_simulation


def generate_bam(
    simulation_params,
    chromosomes,
    include_nonchromosomal,
    simulated_bam_filename,
    genome_fasta_filename,
    simulated_table_filename,
    raw_data_dir,
):
    workflow = pypeliner.workflow.Workflow(default_ctx={'mem':4})

    workflow.setobj(mgd.TempOutputObj('simulation.params'), simulation_params)
    workflow.setobj(mgd.TempOutputObj('chromosomes'), chromosomes)
    workflow.setobj(mgd.TempOutputObj('include_nonchromosomal'), include_nonchromosomal)

    workflow.transform(
        name='create_genome',
        func=destruct.benchmark.destruct_test.create_genome,
        args=(
            mgd.TempInputObj('chromosomes'),
            mgd.TempInputObj('include_nonchromosomal'),
            mgd.OutputFile(genome_fasta_filename),
        ),
    )

    workflow.transform(
        name='create_sim',
        func=destruct.benchmark.create_breakpoint_simulation.create,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.InputFile(genome_fasta_filename),
            mgd.OutputFile(os.path.join(raw_data_dir, 'simulated.fasta')),
            mgd.OutputFile(simulated_table_filename),
            mgd.TempOutputFile('concordant.1.fastq'),
            mgd.TempOutputFile('concordant.2.fastq'),
            mgd.TempOutputFile('discordant.1.fastq'),
            mgd.TempOutputFile('discordant.2.fastq'),
        ),
    )

    workflow.commandline(
        name='cat1',
        args=(
            'cat',
            mgd.TempInputFile('concordant.1.fastq'),
            mgd.TempInputFile('discordant.1.fastq'),
            '>', mgd.OutputFile(os.path.join(raw_data_dir, 'simulated.1.fastq')),
        ),
    )

    workflow.commandline(
        name='cat2',
        args=(
            'cat',
            mgd.TempInputFile('concordant.2.fastq'),
            mgd.TempInputFile('discordant.2.fastq'),
            '>', mgd.OutputFile(os.path.join(raw_data_dir, 'simulated.2.fastq')),
        ),
    )

    workflow.subworkflow(
        name='bwa_align',
        func=destruct.benchmark.align.bwa.workflow.bwa_align_workflow,
        args=(
            mgd.InputFile(genome_fasta_filename),
            mgd.InputFile(os.path.join(raw_data_dir, 'simulated.1.fastq')),
            mgd.InputFile(os.path.join(raw_data_dir, 'simulated.2.fastq')),
            mgd.TempOutputFile('simulated.unsorted.bam'),
        ),
    )

    workflow.transform(
        name='samtools_sort_index',
        func=destruct.benchmark.destruct_test.samtools_sort_index,
        args=(
            mgd.TempInputFile('simulated.unsorted.bam'),
            mgd.OutputFile(simulated_bam_filename),
        ),
    )

    return workflow


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('sim_config',
                           help='Simulation configuration filename')

    argparser.add_argument('ref_data_dir',
                           help='Reference genomes directory')

    argparser.add_argument('raw_data_dir',
                           help='Output directory')

    argparser.add_argument('bam_filename',
                           help='Output bam filename')

    argparser.add_argument('genome_fasta',
                           help='Output genome fasta filename')

    argparser.add_argument('simulated_table',
                           help='Output table of breakpoints in TSV format')

    args = vars(argparser.parse_args())

    try:
        os.makedirs(args['raw_data_dir'])
    except OSError:
        pass

    args['tmpdir'] = os.path.join(args['raw_data_dir'], 'tmp')

    pyp = pypeliner.app.Pypeline(config=args)

    sim_config = yaml.load(open(args['sim_config']))

    ctx = {'mem':4}

    workflow = generate_bam(
        sim_config['simulation'],
        sim_config['reference']['chromosomes'],
        sim_config['reference']['include_nonchromosomal'],
        args['bam_filename'],
        args['genome_fasta'],
        args['simulated_table'],
        args['raw_data_dir'],
    )

    pyp.run(workflow)
