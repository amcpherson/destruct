################################################
# Default configuration for destruct/demix
################################################

import pkg_resources

def get_config(ref_data_dir, user_config):

    ###
    # Gene annotations and reference genome retrieval
    ###

    # Version of ensembl for gene annotations
    ensembl_version                             = '70'

    # Associated genome version used by the ensembl version
    ensembl_genome_version                      = 'GRCh37'

    # Ensemble assemblies to include in the reference genome
    ensembl_assemblies                          = ['chromosome.1', 'chromosome.2', 'chromosome.3', 'chromosome.4', 'chromosome.5', 'chromosome.6', 'chromosome.7', 'chromosome.8', 'chromosome.9', 'chromosome.10', 'chromosome.11', 'chromosome.12', 'chromosome.13', 'chromosome.14', 'chromosome.15', 'chromosome.16', 'chromosome.17', 'chromosome.18', 'chromosome.19', 'chromosome.20', 'chromosome.21', 'chromosome.22', 'chromosome.X', 'chromosome.Y', 'chromosome.MT', 'nonchromosomal']

    # Base chromosomes, used for parallelization and by demix
    chromosomes                                 = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']

    # Ucsc genome version (must match ensembl version!)
    ucsc_genome_version                         = 'hg19'

    # Database of genomic variants genome version (must match ensembl version!)
    dgv_genome_version                          = 'GRCh37_hg19'

    # Database of genomic variants release
    dgv_version                                 = '2013-07-23'

    # Name of the mitochondrial chromosome in ensembl
    mitochondrial_chromosome                    = 'MT'

    ###
    # URLs of gene annotations and reference genome
    ###

    ensembl_assembly_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/fasta/homo_sapiens/dna/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.{0}.fa.gz'

    ensembl_gtf_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/gtf/homo_sapiens/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.gtf.gz'

    dgv_url = 'http://dgv.tcag.ca/dgv/docs/'+dgv_genome_version+'_variants_'+dgv_version+'.txt'

    if ucsc_genome_version == 'hg18':
        rmsk_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/*_rmsk.txt.gz'
    else:
        rmsk_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/'+ucsc_genome_version+'/database/rmsk.txt.gz'

    ###
    # Filenames of gene annotations and reference genome
    ###

    genome_fasta                                = ref_data_dir+'/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.chromosomes.fa'
    genome_fai                                  = genome_fasta+'.fai'
    gtf_filename                                = ref_data_dir+'/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.gtf'
    dgv_filename                                = ref_data_dir+'/dgv.txt'
    repeat_regions                              = ref_data_dir+'/repeats.regions'
    satellite_regions                           = ref_data_dir+'/repeats.satellite.regions'

    # Mapping between ensembl and ucsc chromosome names, hg19 and hg18 provided for you
    chromosome_map                              = pkg_resources.resource_filename('destruct', ucsc_genome_version+'_chr_map.tsv')

    ###
    # Algorithm parameters
    ###

    # Maximum inferred fragment length of a read pair classified as concordant
    bam_max_fragment_length                     = 1000

    # Maximum soft clipped bases before a read is called discordant
    bam_max_soft_clipped                        = 8

    # Maximum fragment length used during realignemnt
    fragment_length_max                         = 800

    # Number of standard deviations to define min and max fragment length more precisely during realignment
    fragment_length_num_stddevs                 = 6.0

    # Minimum alignment probability for filtering discordant read alignments
    alignment_threshold                         = 0.1

    # Prior probability a read is chimeric
    chimeric_prior                              = 0.01

    # Minimum chimeric probability for filtering discordant read alignments
    chimeric_threshold                          = 0.1

    # Minimum read valid probability for filtering discordant read alignments
    readvalid_threshold                         = 0.01

    # Number of concordant reads sampled to calculate valid alignment score distribution
    num_read_samples                            = 100000

    # Minimum width of region covered by discordant reads on either side of the breakpoint
    cluster_coverage_threshold                  = 100

    # Minimum aggregate alignment probability for filtering clusters
    cluster_align_threshold                     = 0.1

    # Minimum aggregate chimeric probability for filtering clusters
    cluster_chimeric_threshold                  = 0.5

    # Minimum aggregate read valid probability for filtering clusters
    cluster_valid_threshold                     = 0.5

    # Minimum discordant read count for filtering clusters
    cluster_readcount_threshold                 = 2

    # Minimum template length aligned to either side of the breakpoint
    template_length_min_threshold               = 40

    # Minimum mate score for possibly concordant alignments
    mate_score_threshold                        = 60

    # Realignment parameters
    match_score                                 = 2
    mismatch_score                              = -3
    gap_score                                   = -4

    # Minimum score for a complex rearrangement cycle
    cycles_scoremax                             = 4

    # Minimum number of breakpoints to visit in search for complex rearrangement cycles
    cycles_visitmax                             = 100000

    # Lambda parameter for complex rearrangement discovery (see nFuse, McPherson et al. 2012)
    cycles_lambda                               = 2000

    ###
    # Demix parameters
    ###

    mappability_length                          = 100
    mappability_filename                        = ref_data_dir+'/'+ucsc_genome_version+'.'+str(mappability_length)+'.bwa.mappability'

    # Thousand genomes dataset
    thousand_genomes_impute_url                 = 'http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz'
    thousand_genomes_dir                        = ref_data_dir+'/ALL_1000G_phase1integrated_v3_impute'
    sample_filename                             = thousand_genomes_dir+'/ALL_1000G_phase1integrated_v3.sample'
    legend_template                             = thousand_genomes_dir+'/ALL_1000G_phase1integrated_v3_chr{0}_impute.legend.gz'
    haplotypes_template                         = thousand_genomes_dir+'/ALL_1000G_phase1integrated_v3_chr{0}_impute.hap.gz'
    genetic_map_template                        = thousand_genomes_dir+'/genetic_map_chr{0}_combined_b37.txt'
    phased_chromosome_x                         = 'X_nonPAR'

    # All snps from thousand genomes
    snp_positions                               = ref_data_dir+'/thousand_genomes_snps.tsv'

    # Heterozygous snp calling
    sequencing_base_call_error                  = 0.01
    het_snp_call_threshold                      = 0.9

    # Shapeit haplotype block resolution
    shapeit_num_samples                         = 100
    shapeit_confidence_threshold                = 0.95

    ###
    # Simulation test parameters
    ###

    reads_per_job                               = 1000000

    ###
    # Parallelization parameters
    ###

    # Number of reads per parallel realignment job
    reads_per_split                             = 1000000

    # Number of clusters per parallel 
    clusters_per_split                          = 1000

    config = locals()
    del config['user_config']

    config.update(user_config)

    return config

