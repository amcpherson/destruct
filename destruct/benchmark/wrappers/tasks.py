import pypeliner.commandline


def concatenate_bcf(in_files, out_file):
    """ Fast concatenation of BCF file using `bcftools`.
    """
    
    cmd = ['bcftools', 'concat', '-a', '-O', 'b', '-o', out_file]
    cmd += [in_files[x] for x in sorted(in_files.keys())]
    
    pypeliner.commandline.execute(*cmd)
                
