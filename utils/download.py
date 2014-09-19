import os
import pypeliner

def wget_gunzip(url, filename):
    temp_filename = filename + '.tmp'
    pypeliner.commandline.execute('wget', url, '-c', '-O', temp_filename + '.gz')
    pypeliner.commandline.execute('gunzip', temp_filename + '.gz')
    os.rename(temp_filename, filename)


def wget(url, filename):
    temp_filename = filename + '.tmp'
    pypeliner.commandline.execute('wget', url, '-c', '-O', temp_filename)
    os.rename(temp_filename, filename)


# REFACTOR: repeated in wrappers.utils
def remove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != 2:
            raise


default_chromosomes = ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT')

ensembl_version = '70'
ensembl_genome_version = 'GRCh37'
default_ensembl_assembly_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/fasta/homo_sapiens/dna/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.{0}.fa.gz'


def download_genome_fasta(genome_fasta,
                          chromosomes=default_chromosomes,
                          include_nonchromosomal=True,
                          assembly_url=default_ensembl_assembly_url):

    ensembl_assemblies = ['chromosome.'+a for a in chromosomes]

    if include_nonchromosomal:
        ensembl_assemblies += ['nonchromosomal']

    temp_filenames = list()

    with open(genome_fasta, 'w') as genome_file:

        for assembly in ensembl_assemblies:

            assembly_fasta = genome_fasta + '.dna.assembly.{0}.fa'.format(assembly)

            if not os.path.exists(assembly_fasta):
                wget_gunzip(assembly_url.format(assembly), assembly_fasta)

            with open(assembly_fasta, 'r') as assembly_file:
                for line in assembly_file:
                    if line[0] == '>':
                        line = line.split()[0] + '\n'
                    genome_file.write(line)

            temp_filenames.append(assembly_fasta)

    for temp_filename in temp_filenames:
        remove(temp_filename)


