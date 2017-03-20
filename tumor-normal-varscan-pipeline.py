'''
Pipeline to analyse tumor-normal WES paired data from fastq.
'''

from collections import defaultdict
import argparse
import ConfigParser
import gzip
import itertools
import os
import jinja2
import re
import shutil
import subprocess
import sys


def align_reads(fastq_tsv, genome, config, dir_map):
    '''
    Aligns reads from fastqs.
    '''
    bwa = config.get('Binaries', 'bwa')
    samtools = config.get('Binaries', 'samtools')
    ref_genome = config.get(genome, 'ref_genome')
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
        subprocess.call(cmd, shell=True)
        bams[sample_name].append(bam)
    return bams


def call_variants(pileups, contrasts, config, dir_map):
    '''
    From pileups and contrasts, varscan2 somatic for variant calling.
    '''
    java = config.get('Binaries', 'java')
    varscan = config.get('Binaries', 'varscan')
    bcftools = config.get('Binaries', 'bcftools')
    vcf_map = {}
    for normal_samplename, tumor_samplename in contrasts:
        normal = pileups[normal_samplename]
        tumor = pileups[tumor_samplename]
        samplename = "_vs_".join([normal_samplename, tumor_samplename])
        snps = '/'.join([dir_map["vcfdir"],
                         samplename + config.get('Suffix', 'snps')])
        indels = '/'.join([dir_map["vcfdir"],
                           samplename + config.get('Suffix', 'indels')])
        vcf = '/'.join([dir_map["vcfdir"],
                        samplename + config.get('Suffix', 'vcf')])
        tmp = '/'.join([dir_map["vcfdir"],
                        "tmp.vcf"])
        cmd1 = ' '.join([java, "-jar", varscan, "somatic", normal, tumor,
                        samplename, "--p-value 0.99 --output-vcf 1"])
        cmd2 = ' '.join([bcftools, "concat", snps, indels, ">", tmp])
        # need to change
        normal_name = normal_samplename + ':' + samplename
        tumor_name = tumor_samplename + ':' + samplename
        rename_header_samplenames(tmp, vcf, normal_name, tumor_name)
        subprocess.call(cmd1, shell=True)
        subprocess.call(cmd2, shell=True)
        vcf_map[tumor_samplename] = snps
    return vcf_map


def rename_header_samplenames(vcf, new_vcf, normal_name, tumor_name):
    '''
    Renames the samplename headers of the VCF.
    '''
    is_renamed = False
    with open(vcf, 'rU') as vcf_handle, open(new_vcf, 'w') as new_handle:
        for line in vcf_handle:
            if not is_renamed and re.match("#CHROM", line):
                arow = line.strip('\n').split
                arow[9] = normal_name
                arow[10] = tumor_name
                new_handle.write('\t'.join(arow) + '\n')
            else:
                new_handle.write(line)


def create_pileups(bam_map, genome, config, dir_map):
    '''
    Creates pileups from bam files with samtools mpileup.
    '''
    samtools = config.get('Binaries', 'samtools')
    ref_genome = config.get(genome, 'ref_genome')
    pileup_map = {}
    for samplename, bam in bam_map.items():
        pileup = '/'.join([dir_map["pileupdir"],
                           samplename + config.get('Suffix', 'pileup')])
        cmd = ' '.join([samtools, "mpileup -B -f", ref_genome, bam,
                        ">", pileup])
        subprocess.call(cmd, shell=True)
        pileup_map[samplename] = pileup
    return pileup_map


def get_rg_values(sample_name, fastq):
    '''
    Parses gzipped fastq for RG values and returns RG str for BWA.
    '''
    with gzip.open(fastq, 'rU') as handle:
        line = handle.readline()
        arow = line.strip('\n').split()
        info = arow[0].split(':')[1:]
        instrument_id = info[0]
        run_id = info[1]
        flowcell_id = info[2]
        flowcell_lane = info[3]
        index_seq = arow[1].split(':')[3]
        rgid = '.'.join([sample_name, flowcell_id, flowcell_lane])
    rglb = '.'.join([sample_name, run_id])
    rgpu = '.'.join([instrument_id,
                     flowcell_lane,
                     index_seq])
    rgsm = sample_name
    rgcn = "DFCI-CCCB"
    rgpl = "ILLUMINA"
    rg_vals = "@RG\\tID:" + rgid + "\\tPL:" + rgpl + "\\tLB:" + \
              rglb + "\\tSM:" + rgsm + "\\tCN:" + rgcn + "\\tPU:" + rgpu
    return rg_vals


def get_sample_pairs(contrasts):
    '''
    From two column TSV, get matched pairs.
    '''
    pairs = []
    with open(contrasts, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            pairs.append(arow)
    return pairs


def map_fastq_from_tsv(fastq_tsv):
    '''
    Parses tsv into map of sample name and lane specific fastq.
    '''
    fastq_map = defaultdict(list)
    with open(fastq_tsv, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            sample = arow[0]
            fastq = arow[1]
            fastq_map[sample].append(fastq)
    return fastq_map


def map_fastq_pairs(fastqs):
    '''
    Gets fastqs from fastq dir and returns dict (k=first, v=second).
    '''
    first_read_fastqs = [f for f in fastqs if re.search("_R1_", f)]
    second_read_fastqs = [f for f in fastqs if re.search("_R2_", f)]
    basename_second_map = {f.replace("_R2_", "") : f 
                           for f in second_read_fastqs}
    mapped_pairs = {f : basename_second_map[f.replace("_R1_", "")]
                    for f in first_read_fastqs
                    if f.replace("_R1_", "") in basename_second_map}
    first_only = {f : "" for f in first_read_fastqs
                  if f.replace("_R1_", "") not in basename_second_map}
    second_only = {f : "" for f in second_read_fastqs
                   if f.replace("_R2_", "" not in mapped_pairs)}
    mapped_pairs.update(first_only)
    mapped_pairs.update(second_only)
    return mapped_pairs


def merge_bams(indv_bams_map, config, dir_map,dry=False):
    '''
    Merges the individual bams into a sample bam.
    '''
    merged_bams = {}
    samtools = config.get('Binaries', 'samtools')
    for samplename, bams in indv_bams_map.items():
        input_bams = ' '.join(bams)
        output_bam = '/'.join([dir_map["bamdir"], 
                               samplename + config.get('Suffix', 'bam')])
        cmd1 = ' '.join([samtools, "merge", output_bam, input_bams])
        cmd2 = ' '.join([samtools, "index", output_bam])
        if dry:
            sys.stdout.write(cmd1 + '\n' + cmd2 + '\n')
        else:
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
        merged_bams[samplename] = output_bam
    return merged_bams


def merge_vcf(vcf_map, config, dir_map, dry=False):
    '''
    Merges the VCFs into one for putting into VEP and gemini.
    '''
    bcftools = config.get('Binaries', 'bcftools')
    merged_vcf = '/'.join([dir_map["vcfdir"], "all_merged.vcf"])
    vcfs = ' '.join(vcf_map.values())
    cmd = ' '.join([$bcftools])
    return merged_vcf


def realign_indels(bam_map, genome, config, dir_map):
    '''
    Realigns indels with GATK.
    '''
    java = config.get('Binaries', 'java')
    gatk = config.get('Binaries', 'gatk')
    ref_genome = config.get(genome, 'ref_genome')
    realn_bams = {}
    for samplename, bam in bam_map.items():
        realn_intervals = dir_map["bamdir"] + '/' + samplename + \
                          config.get('Suffix', 'indelrealnintervals')
        realn_bam = dir_map["bamdir"] + '/' + samplename + \
                    config.get('Suffix', 'indelrealn')
        cmd1 = ' '.join([java, "-jar", gatk, "-T RealignerTargetCreator",
                         "-R", ref_genome, "-I", bam, "-o", realn_intervals])
        cmd2 = ' '.join([java, "-jar", gatk, "-T IndelRealigner",
                         "-R", ref_genome, "-I", bam,
                         "-targetIntervals", realn_intervals, "-o", realn_bam])
        subprocess.call(cmd1, shell=True)
        subprocess.call(cmd2, shell=True)
        realn_bams[samplename] = realn_bam
    return realn_bams


def setup_dir(cur_dir, out_dir_name, dry=False):
    '''
    Sets up output directory in project directory.
    '''
    if cur_dir == '.':
        cur_dir = os.getcwd()
    if not os.path.isdir(cur_dir):
        sys.stderr.write("Error: project directory path does not exist.\n")
        sys.exit()
    out_dir = '/'.join([cur_dir, out_dir_name])
    bam_dir = '/'.join([out_dir, "bam_files"])
    indbam_dir = '/'.join([bam_dir, "individual_bam_files"])
    pileup_dir = '/'.join([out_dir, "pileups"])
    vcf_dir = '/'.join([out_dir, "vcfs"])
    report_dir = '/'.join([out_dir, "report_html"])
    coverage_dir = '/'.join([out_dir, "coverage"])
    if not dry:
        try:
            for folder in [out_dir, bam_dir, indbam_dir, pileup_dir,
                           coverage_dir, vcf_dir, report_dir]:
                os.makedirs(folder)
        except OSError as err:
            sys.stderr.write("%s\n" % err)
            sys.stderr.write("Error: %s directory already exists.\n" % folder)
            sys.exit()
    return {"bamdir": bam_dir,
            "outdir": out_dir,
            "projdir": cur_dir,
            "indbamdir": indbam_dir,
            "pileupdir": pileup_dir,
            "vcfdir": vcf_dir,
            "coveragedir": coverage_dir,
            "reportdir": report_dir}


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
    # Set up globally used maps.
    dir_map = setup_dir(args.projectdir, args.outdir)
    config = ConfigParser.ConfigParser()
    script_dir = os.path.dirname(os.path.realpath(__file__))
    config.read(os.path.join(script_dir, "config"))
    # Processing
    indv_bams_maps = align_reads(args.fastq_tsv, args.genome, config, dir_map)
    merged_bams_map = merge_bams(indv_bams_maps, config, dir_map)
    realn_bams_map = realign_indels(merged_bams_map, args.genome, config,
                                    dir_map)
    pileups_map = create_pileups(realn_bams_map, args.genome, config, dir_map)
    sample_pair_maps = get_sample_pairs(args.contrasts)
    vcf_map = call_variants(pileups_map, sample_pair_maps,
                                               config, dir_map)
    merge_vcf(vcf_map, sample_pa)


if __name__ == "__main__":
    main()
