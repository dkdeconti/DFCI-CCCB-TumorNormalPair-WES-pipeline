#!/usr/bin/python
'''
Main entry for plexwell amplicon analysis.
'''

from collections import defaultdict
import argparse
import ConfigParser
import gzip
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import jinja2
import re
import shutil
import subprocess
import sys


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


def call_variants(pileups_map, config, dir_map, dry=False):
    '''
    Pipes samtools mpileup of bam to varscan.
    '''
    java = config.get('Binaries', 'java')
    varscan = config.get('Binaries', 'varscan')
    snp_map = {}
    indel_map = {}
    for samplename, pileup in pileups_map.items():
        snps = '/'.join([dir_map["vcfdir"],
                         samplename + config.get('Suffix', 'snps')])
        indels = '/'.join([dir_map["vcfdir"],
                           samplename + config.get('Suffix', 'indels')])
        cmd1 = ' '.join([java, "-jar", varscan, "mpileup2snp", pileup,
                         "--p-value 99e-02 --output-vcf 1", ">", snps])
        cmd2 = ' '.join([java, "-jar", varscan, "mpileup2indel", pileup,
                         "--p-value 99e-02 --output-vcf 1", ">", indels])
        if dry:
            sys.stdout.write(cmd1 + '\n' + cmd2 + '\n')
        else:
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
        snp_map[samplename] = snps
        indel_map[samplename] = indels
    return snp_map, indel_map


def cluster_filter(variant_map, cluster_map):
    '''
    Filters variants by intersection with clusters.
    '''
    filtered_variant_map = {}
    for samplename, variants in variant_map.items():
        clusters = cluster_map[samplename]
        filtered_variants = [v for i, v in enumerate(variants)
                             if i == 0 or has_intersect(v, clusters)]
        filtered_variant_map[samplename] = filtered_variants
    return filtered_variant_map


def cluster_regions(bams_map, min_mean_depth, config, dir_map, dry=False):
    '''
    Determines specific amplicon by min_mean_coverage.
    '''
    bedtools = config.get('Binaries', 'bedtools')
    bed_map = {}
    cluster_map = {}
    for samplename, bam in bams_map.items():
        bed = '/'.join([dir_map["coveragedir"],
                        samplename + config.get('Suffix', 'coveragebed')])
        cmd = ' '.join([bedtools, "merge -i", bam, "| head -n -1 |",
                        bedtools, "coverage -a - -b", bam, ">", bed])
        if dry:
            sys.stdout.write(cmd + '\n')
        else:
            subprocess.call(cmd, shell=True)
        bed_map[samplename] = bed
        cluster_map[samplename] = filter_clusters(bed, min_mean_depth)
    return bed_map, cluster_map


def create_report(snps, indels, plots, dir_map, dry=False):
    '''
    Injects data into html template with jinja2.
    '''
    this_dir = os.path.dirname(os.path.realpath(__file__))
    lib_dir = os.path.join(this_dir, 'lib')
    report_dir = dir_map["reportdir"]
    lib_destination = os.path.join(report_dir, 'lib')
    report = '/'.join([report_dir, "html_report.html"])
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(this_dir))
    # can't seem to find template.html
    template = env.get_template("template.html")
    samples = {samplename : {"header" : snps[samplename][0],
                             "snps" : snps[samplename][1:],
                             "indels" : indels[samplename][1:],
                             "plots": plots[samplename],
                             "samplename": samplename.split('/')[-1]}
               for samplename in snps.keys()}
    context = {"samples": samples}
    if not dry:
        with open(report, 'w') as outfile:
            outfile.write(template.render(context))
        shutil.copytree(lib_dir, lib_destination)


def filter_clusters(bed, min_depth):
    '''
    Parses BED file to filter for high coverage regions.
    '''
    #clusters = defaultdict(list)
    clusters = {} # for some reason defaultdict is buggy
    with open(bed, 'rU') as bedfile:
        for line in bedfile:
            arow = line.strip('\n').split('\t')
            chrom = arow[0]
            try:
                begin = int(arow[1])
                end = int(arow[2])
                total_depth = int(arow[3])
                bases_covered = int(arow[5])
            except ValueError as err:
                sys.stderr.write(str(err) + "\nError in coverage bed file.\n")
                sys.exit()
            if total_depth >= min_depth:
                #clusters[chrom].append((begin, end))
                if chrom in clusters:
                    clusters[chrom].append((begin, end))
                else:
                    clusters[chrom] = [(begin, end)]
    return clusters


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


def get_matched_cluster(pos, clusters):
    '''
    Produces a key based on which cluster pos resides.
    '''
    vk = pos[0]
    vpos = int(pos[1])
    for cluster in clusters[vk]:
        begin = int(cluster[0])
        end = int(cluster[1])
        #print "DEBUG VALUES", vk, vpos, begin, end
        #print "DEBUG TYPES", type(vk), type(vpos), type(begin), type(end)
        if begin <= vpos <= end:
            return ("%s:%i-%i" % (vk, begin, end), begin)


def has_intersect(variant, clusters):
    '''
    Checks for any intersects between variant and any cluster.
    '''
    intersected = False
    vk = variant["chrom"]
    try:
        vpos = int(variant["pos"])
    except ValueError as err:
        sys.stderr.write(str(err) + "\n")
        sys.exit()
    return_bool = False
    if vk in clusters:
        for cluster in clusters[vk]:
            begin = int(cluster[0])
            end = int(cluster[1])
            if begin <= vpos <= end:
                return_bool = True
                break
    return return_bool
    # elegance not working...
    #return vk in clusters or (len(clusters[vk]) > 0 
    # and any(b[0] <= vpos <= b[1] for b in clusters[vk]))


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


def parse_pileup_for_qc_stats(pileup_map, cluster_map, dir_map, dry=False):
    '''
    Parses pileup to output pileup stats.
    '''
    stats = {}
    for samplename, pileup in pileup_map.items():
        clusters = cluster_map[samplename]
        stats[samplename] = {}
        with open(pileup, 'rU') as pileup_file:
            for line in pileup_file:
                arow = line.strip('\n').split('\t')
                chrom = arow[0]
                pos = int(arow[1])
                if has_intersect({"chrom": chrom, "pos": pos}, clusters):
                    cluster_key, start = get_matched_cluster((chrom, pos),
                                                             clusters)
                else:
                    continue
                if cluster_key not in stats[samplename]:
                    stats[samplename][cluster_key] = [[], [], [], [], []]
                adj_pos = pos - start
                depth = int(arow[3])
                qual = [ord(c) - 33 for c in arow[5]]
                mq = np.median(qual)
                uq = np.percentile(qual, 90)
                lq = np.percentile(qual, 10)
                for i, v in enumerate([adj_pos, depth, mq, uq, lq]):
                    stats[samplename][cluster_key][i].append(v)
    return stats


def parse_vcf(vcf, is_indel=False, dry=False):
    '''
    Parses VCFs to a map of select features.
    '''
    header = {"chrom": "chrom",
              "pos": "position",
              "ref": "ref allele",
              "alt": "alt allele",
              "freq": "alt freq",
              "pval": "p value",
              "total_depth": "total depth",
              "ref_depth": "ref depth",
              "allele_depth": "alt depth"}
    variants = [header]
    with open(vcf, 'rU') as handle:
        for line in handle:
            if line[0] == "#":
                continue # skip header
            vrow = line.strip('\n').split('\t')
            chrom = vrow[0]
            pos = vrow[1]
            if is_indel:
                pos = str(int(vrow[1]) + 1)
            ref = vrow[3]
            alt = vrow[4]
            v_info = vrow[9].split(':')
            freq = v_info[6]
            pval = v_info[7]
            total_depth = v_info[3]
            ref_depth = v_info[4]
            allele_depth = v_info[5]
            variant = {"chrom": chrom,
                       "pos": pos,
                       "ref": ref,
                       "alt": alt,
                       "freq": freq,
                       "pval": pval,
                       "total_depth": total_depth,
                       "ref_depth": ref_depth,
                       "allele_depth": allele_depth}
            variants.append(variant)
    return variants


def pileup(bam_map, config, dir_map, dry=False):
    '''
    Creates pileups from bam files with samtools mpileup.
    '''
    samtools = config.get('Binaries', 'samtools')
    ref_genome = config.get('References', 'ref_genome')
    pileup_map = {}
    for samplename, bam in bam_map.items():
        pileup = '/'.join([dir_map["pileupdir"],
                           samplename + config.get('Suffix', 'pileup')])
        cmd = ' '.join([samtools, "mpileup -B -f", ref_genome, bam,
                        ">", pileup])
        if dry:
            sys.stdout.write(cmd + "\n")
        else:
            subprocess.call(cmd, shell=True)
        pileup_map[samplename] = pileup
    return pileup_map


def plot_qc(stats_map, config, dir_map, dry=False):
    '''
    Plots depth and base quality across regions.
    '''
    plot_map = defaultdict(list)
    for samplename, stats in stats_map.items():
        for cluster_key, stat_matrix in stats.items():
            filename = samplename + "." + cluster_key + \
                       config.get('Suffix', 'qcstats')
            print "DEBUG FILENAME", filename
            plot_file = '/'.join([dir_map["reportdir"], filename])
            pos = stat_matrix[0]
            depth = stat_matrix[1]
            mq = stat_matrix[2]
            uq = stat_matrix[3]
            lq = stat_matrix[4]
            # Start figure making
            fig, host = plt.subplots(figsize=(16, 8))
            par1 = host.twinx()
            p1 = host.plot(pos, mq, "b-", label="Base Quality")
            p2 = par1.plot(pos, depth, "r-", label="Depth")
            host.set_xlabel("Position along %s" % cluster_key)
            host.set_ylabel("Base quality score")
            par1.set_ylabel("Depth")
            host.yaxis.label.set_color("b")
            par1.yaxis.label.set_color("r")
            host.tick_params(axis="y", colors="b")
            par1.tick_params(axis="y", colors="r")
            lines = [p1, p2]
            plt.savefig(plot_file)
            plt.close()
            plot_map[samplename].append(plot_file)
    return plot_map


def realign_indels(bam_map, config, dir_map, dry=False):
    '''
    Realigns indels with GATK.
    '''
    java = config.get('Binaries', 'java')
    gatk = config.get('Binaries', 'gatk')
    ref_genome = config.get('References', 'ref_genome')
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
        if dry:
            sys.stdout.write(cmd1 + '\n' + cmd2 + '\n')
        else:
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
    Parses CLI args and central dispatch of functions.
    '''
    # Arg parsing
    parser = argparse.ArgumentParser(description="Plexwell amplicon analysis")
    parser.add_argument("fastqtsv", metavar="FastqTSV",
                        help="lane specific fastq TSV")
    parser.add_argument("-d", "--projectdir",
                        metavar="PROJECT_DIR",
                        help="project directory; default: .")
    parser.add_argument("-o", "--outdir",
                        metavar="OUTPUT_DIR",
                        help="output dir; relative to -d; default: plexout")
    parser.add_argument("-g", "--genome",
                        help="genome [hg19,]; default: hg19")
    parser.add_argument("-m", "--mindepth",
                        type=int,
                        metavar="MIN_MEAN_DEPTH",
                        help="min mean coverage for clusters; default 200")
    parser.add_argument("--dry",
                        action="store_true",
                        help="Creates dirs, but no file creation")
    parser.set_defaults(projectdir=".", outdir="plexout", genome="hg19",
                        mindepth=5000)
    args = parser.parse_args()
    # Set up globally used maps.
    bin_map = {"samtools": "~/bin/samtools",
            "bwa": "~/bin/bwa",
            "java": "~/bin/java",
            "gatk": "~/bin/gatk.jar",
            "varscan": "~/bin/varscan.jar",
            "bedtools": "~/bin/bedtools"}
    suffix_map = {"indelrealn": ".indelrealn.bam",
                  "indelrealnintervals": ".indelrealn.intervals",
                  "coveragebed": ".bed",
                  "bam": ".bam",
                  "snps": ".snps.vcf",
                  "indels": ".indels.vcf",
                  "pileup": ".pileup",
                  "qcstats": ".qcstats.png"}
    dir_map = setup_dir(args.projectdir, args.outdir, dry=args.dry)
    config = ConfigParser.ConfigParser()
    script_dir = os.path.dirname(os.path.realpath(__file__))
    config.read(os.path.join(script_dir, "config"))
    # Processing
    indv_bams_map = align_reads(args.fastqtsv, config, dir_map, dry=args.dry)
    merged_bams_map = merge_bams(indv_bams_map, config, dir_map, dry=args.dry)
    realn_bams_map = realign_indels(merged_bams_map, config, dir_map,
                                    dry=args.dry)
    _, cluster_map = cluster_regions(realn_bams_map, int(args.mindepth),
                                     config, dir_map, dry=args.dry)
    pileups_map = pileup(realn_bams_map, config, dir_map, dry=args.dry)
    stats_map = parse_pileup_for_qc_stats(pileups_map, cluster_map, dir_map,
                                          dry=args.dry)
    plots_map = plot_qc(stats_map, config, dir_map, dry=args.dry)
    snp_vcf_map, indel_vcf_map = call_variants(pileups_map, config, dir_map,
                                               dry=args.dry)
    snp_map = {samplename : parse_vcf(vcf, dry=args.dry)
               for samplename, vcf in snp_vcf_map.items()}
    indel_map = {samplename : parse_vcf(vcf, is_indel=True, dry=args.dry)
                 for samplename, vcf in indel_vcf_map.items()}
    f_snp_map, f_indel_map = [cluster_filter(m, cluster_map)
                              for m in [snp_map, indel_map]]
    # Create html report.
    create_report(f_snp_map, f_indel_map, plots_map, dir_map)


if __name__ == "__main__":
    main()
