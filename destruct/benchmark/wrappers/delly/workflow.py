import os
import shutil

import pypeliner
import pypeliner.managed as mgd

import destruct.benchmark.wrappers.tasks
import destruct.benchmark.wrappers.delly.tasks
import destruct.benchmark.wrappers.utils


def create_delly_wrapper_workflow(bam_filenames, output_filename, raw_data_dir, control_id=None, ref_genome_fasta_file=None, delly_excl_chrom=None):
    bams = list()
    for lib_id, bam_filename in bam_filenames.items():
        bams += [destruct.benchmark.wrappers.utils.symlink(bam_filename, link_name='{0}.bam'.format(lib_id), link_directory=raw_data_dir)]
        destruct.benchmark.wrappers.utils.symlink(bam_filename+'.bai', link_name='{0}.bam.bai'.format(lib_id), link_directory=raw_data_dir)

    workflow = pypeliner.workflow.Workflow()
    
    workflow.transform(
        name='get_sv_types',
        func=destruct.benchmark.wrappers.delly.tasks.get_sv_types,
        ret=pypeliner.managed.OutputChunks('sv_type'),
        args=(
            mgd.InputFile(ref_genome_fasta_file),
        ),
    )

    workflow.transform(
        name='delly_call',
        axes=('sv_type',),
        ctx={'mem': 64, 'num_retry': 2, 'mem_retry_factor': 2},
        func=destruct.benchmark.wrappers.delly.tasks.run_delly_call,
        args=(
            mgd.Instance('sv_type'),
            delly_excl_chrom,
            ref_genome_fasta_file,
            [mgd.InputFile(bam) for bam in bams],
            mgd.TempOutputFile('out.bcf', 'sv_type'),
        ),
    )

    if control_id is None:
        concat_input = mgd.TempInputFile('out.bcf', 'sv_type')

    else:
        workflow.transform(
            name='delly_filter_somatic',
            axes=('sv_type',),
            ctx={'mem': 4, 'num_retry': 2, 'mem_retry_factor': 2},
            func=destruct.benchmark.wrappers.delly.tasks.run_delly_filter,
            args=(
                mgd.Instance('sv_type'),
                bam_filenames.keys(),
                control_id, 
                mgd.TempSpace('samples.tsv'),
                ref_genome_fasta_file,
                mgd.TempInputFile('out.bcf', 'sv_type'),
                mgd.TempOutputFile('somatic.bcf', 'sv_type'),
            ),
        )

        concat_input = mgd.TempInputFile('somatic.bcf', 'sv_type')

    workflow.transform(
        name='concatenate_vcf',
        func=destruct.benchmark.wrappers.tasks.concatenate_bcf,
        ctx={'mem': 4, 'num_retry': 2, 'mem_retry_factor': 2},
        args=(
            concat_input,
            mgd.TempOutputFile('somatic.bcf'),
        ),
    )

    workflow.transform(
        name='convert_vcf',
        func=destruct.benchmark.wrappers.delly.tasks.convert_vcf,
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        args=(
            mgd.TempInputFile('somatic.bcf'),
            mgd.OutputFile(output_filename),
        ),
        kwargs={
            'control_id': control_id,
        }
    )

    return workflow


def setup_delly(test_config, **kwargs):
    ensembl_version = kwargs.get('ensembl_version', '70')
    ensembl_genome_version = kwargs.get('ensembl_genome_version', 'GRCh37')

    primary_assembly_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/fasta/homo_sapiens/dna/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.primary_assembly.fa.gz'
    chromosome_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/fasta/homo_sapiens/dna/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.chromosome.{0}.fa.gz'
    hg19_excl_url = 'http://raw.githubusercontent.com/tobiasrausch/delly/master/excludeTemplates/human.hg19.excl.tsv'

    if ensembl_genome_version != 'GRCh37':
        raise Exception('Only hg19 supported')
    
    genome_fasta = kwargs['ref_genome_fasta_file']
    delly_excl_chrom = kwargs['delly_excl_chrom']

    destruct.benchmark.wrappers.utils.makedirs(os.path.dirname(genome_fasta))
    destruct.benchmark.wrappers.utils.makedirs(os.path.dirname(delly_excl_chrom))

    if test_config.get('chromosomes', None) is None:
        destruct.benchmark.wrappers.utils.wget_file_gunzip(primary_assembly_url, genome_fasta)

    else:
        with open(genome_fasta, 'w') as genome_file:
            for chromosome in test_config['chromosomes']:
                chromosome_filename = './chromosome_{0}.fa'.format(chromosome)

                destruct.benchmark.wrappers.utils.wget_file_gunzip(chromosome_url.format(chromosome), chromosome_filename)

                with open(chromosome_filename, 'r') as chromosome_file:
                    shutil.copyfileobj(chromosome_file, genome_file)

                os.remove(chromosome_filename)

    destruct.benchmark.wrappers.utils.wget_file(hg19_excl_url, delly_excl_chrom)
