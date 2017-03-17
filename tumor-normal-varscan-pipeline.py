'''

'''

def align_reads(fastq_tsv, config, dir_map, dry=False):
    '''
    Aligns reads from fastqs.
    '''
    bwa = config.get('Binaries', 'bwa')
    samtools = config.get('Binaries', 'samtools')
    ref_genome = config.get('References', 'ref_genome')
    fastq_map = map_fastq_from_tsv(fastq_tsv)
    inverted_fastq_map = dict((v, k) for k in fastq_map for v in fastq_map[k])
    fastqs = itertools.chain.from_iterable(fastq_map.values())
    fastq_pairs = map_fastq_pairs(fastqs)
    bams = defaultdict(list)
    for first, second in fastq_pairs.items():
        sample_name = inverted_fastq_map[first]
        bam_basename = re.sub(r'_R[1-2]_[0-9]+\.fastq.gz',
                              '', first).split('/')[-1]
        bam = dir_map["indbamdir"] + '/' + bam_basename + \
              config.get('Suffix', 'bam')
        rg_vals = get_rg_values(sample_name, first)
        cmd = ' '.join([bwa, 'mem -R \"%s\"' % rg_vals,
                        ref_genome, first, second, '|',
                        samtools, "view -bht", ref_genome, "|",
                        samtools, "sort", ">", bam])
        if dry:
            sys.stdout.write(cmd + '\n')
        else:
            subprocess.call(cmd, shell=True)
        bams[sample_name].append(bam)
    return bams


def main():
    '''
    Parses CLI args and setup.
    '''
    parser = argparse.ArgumentParser(description="VarScan Paired WES analysis")
    parser.add_argument("fastqtsv", metavar="FastqTSV",
                        help="lane specific fastq TSV")
    parser.add_argument("contrasts", metavar="ContrastFile",
                        help="Two column TSV of pairs")
    parser.add_argument("-d", "--projectdir",
                        metavar="PROJECT_DIR",
                        help="project directory; default: .")
    parser.add_argument("-o", "--outdir",
                        metavar="OUTPUT_DIR",
                        help="output dir; relative to -d; default: paired_WES")
    parser.add_argument("-g", "--genome",
                        help="genome [hg19,]; default: hg19")
    parser.set_defaults(projectdir=".", outdir="paired_WES", genome="hg19")
    args = parser.parse_args()
    
    pass


if __name__ == "__main__":
    main()
