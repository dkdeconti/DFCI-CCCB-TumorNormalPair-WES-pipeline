'''
Converts bedtools coverage -hist bed file to GATK DepthOfCoverage output.
'''
import argparse
import re


def make_headers():
    '''
    Makes header for sample.coverage.
    '''
    hist_header =  '\t'.join(["Source_of_reads"] +
                             ["from_%i_to_%i)" % (i, i+1) for i in range(0, 500)] +
                             ["from_500_to_inf"])
    summary_header = '\t'.join(["sample_id", "total", "mean",
                                "granular_third_quartile",
                                "granular_median",
                                "grandular_first_quartile",
                                "%_bases_above_15"])
    return summary_header, hist_header


def parse_bed(bedfile):
    '''
    Parses bedfile into histogram dict and header.
    '''
    total_depth = 0
    total_bases = 0
    hist = {i : 0 for i in range(0, 501)}
    with open(bedfile, 'rU') as handle:
        for line in handle:
            if re.match("all", line):
                arow = line.strip('\n').split('\t')
                depth = int(arow[1])
                num_bases = int(arow[2])
                total_depth += depth
                total_bases += num_bases
                if depth > 500:
                    hist[500] += num_bases
                else:
                    hist[depth] += num_bases
    mean = total_depth/float(total_bases)
    summ_out = [str(total_depth), str(mean), "Na", "Na", "Na", "Na"]
    hist_out = [str(hist[k]) for k in sorted(hist.keys())]
    return summ_out, hist_out


def main():
    '''
    Parses CLI args and setup.
    '''
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("bed", metavar="BED", nargs='+',
                        help="bedtools coverage -hist BED")
    args = parser.parse_args()
    summary_header, hist_header = make_headers()
    with open('coverage.sample_statistics', 'w') as cov_stats, \
         open('coverage.sample_summary', 'w') as cov_summ:
        cov_stats.write(hist_header + '\n')
        cov_summ.write(summary_header + '\n')
        for bed in args.bed:
            sample_name = bed.split('.')[0]
            summ_out, hist_out = parse_bed(bed)
            cov_stats.write('\t'.join([sample_name] + hist_out) + '\n')
            cov_summ.write('\t'.join([sample_name] + summ_out) + '\n')


if __name__ == "__main__":
    main()
