#!/usr/bin/python
# Since I have to frequently grab this table, output a table of MAPD values for a particular VCF or set of 
# VCF files.
#
# 3/28/2016 - D Sims
###########################################################################################################
import sys
import os
import re
import subprocess
import datetime
import argparse
import subprocess
import math
from pprint import pprint as pp
from distutils.version import LooseVersion

version = '3.5.2_072617'

def get_args():
    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 100, width=200),
        description='''
        Input one or more VCF files and generate a table of some important metrics. Probably going to be removed
        or replaced with better script at some point soon.
        ''',
        version = '%(prog)s  - ' + version,
    )
    parser.add_argument('vcf', metavar='<VCF(s)>', nargs='+', help='VCF file(s) to process')
    parser.add_argument('-d','--dna_only', action='store_true', 
            help='Data comes from DNA only specimens and no RNA data to report. Essentially reporting MAPD only.')
    parser.add_argument('-o', '--output', metavar='<outfile>', help='Custom output file (DEFAULT: %(default)s)')
    args = parser.parse_args()
    return args

def read_vcf(vcf_file):
    oca_v3_version = '2.3'
    expr_sum = 0
    fetched_data = {}

    try:
        with open(vcf_file) as fh:
            for line in fh:
                if not line.startswith('#') and 'SVTYPE=ExprControl' not in line:
                    continue
                elif line.startswith('##mapd='):
                    mapd = get_value(line.rstrip())
                    if float(mapd) > 0.5:
                        mapd = flag_val(mapd)
                    fetched_data['MAPD'] = mapd

                elif line.startswith('##TotalMappedFusionPanelReads='):
                    rna_reads = get_value(line.rstrip())
                    if int(rna_reads) < 100000:
                        rna_reads = flag_val(rna_reads)
                    fetched_data['RNA_Reads'] = rna_reads

                elif line.startswith('##fileDate'):
                    date = datetime.datetime.strptime(get_value(line.rstrip()),"%Y%M%d")
                    formatted_date = date.strftime("%Y-%M-%d")
                    fetched_data['Date'] = (formatted_date)

                elif line.startswith('##OncomineVariantAnnotationToolVersion'):
                    ovat_version = get_value(line.rstrip())
                    if LooseVersion(ovat_version) > LooseVersion(oca_v3_version):
                        p1,p2 = get_rna_pool_info(vcf_file)
                        if int(p1) < 100000:
                            p1 = flag_val(p1)
                        if int(p2) < 100000:
                            p2 = flag_val(p2)
                        fetched_data['Pool1'] = p1
                        fetched_data['Pool2'] = p2

                elif re.search('SVTYPE=ExprControl',line):
                    read_count = re.search('READ_COUNT=(\d+).*',line).group(1)
                    expr_sum += int(read_count)

            if expr_sum < 20000:
                expr_sum = flag_val(str(expr_sum))
            fetched_data['Expr_Sum'] = expr_sum

    except IOError as e:
        sys.stderr.write('ERROR: Can not open file {}: {}!\n'.format(vcf_file,e))
        sys.exit()
    return fetched_data

def flag_val(val):
    return '*' + val + '*'

def get_rna_pool_info(vcf):
    p = subprocess.Popen(['match_rna_qc.pl', vcf], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (data,err) = p.communicate()
    ret_res = data.split('\n')
    results = dict(zip(ret_res[0].split(','),ret_res[1].split(',')))
    p1_tot = str(round(float(results['pool1_total'])))
    p2_tot = str(round(float(results['pool2_total'])))
    return p1_tot.split('.')[0], p2_tot.split('.')[0]

def get_value(line):
    return line.split('=')[1]

def get_name_from_vcf(vcf):
    with open(vcf) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                return line.split('\t')[-1].rstrip('\n')

def col_size(data):
    col_width = 0
    for i in data:
        if len(i) > col_width:
            col_width = len(i)
    return col_width + 4

def print_data(results,outfile,dna_only):
    # Figure out if we have two different versions of analysis, and if so bail out to make easier.
    l = [len(v) for k,v in results.iteritems()]
    if len(set(l)) > 1:
        print('Mixed version VCFs detected!  We can not process two different versions together!  '
              'Please run separately and cat the data later.')
        sys.exit(1)

    # Check to make sure we didn't really mean to use DNA only.
    if not dna_only:
        for sample in results:
            try:
                tmp_reads = results[sample]['RNA_Reads']
            except KeyError:
                print('ERROR: DNA only specimen(s) detected.  You must runs these using the "--dna_only" option!')
                sys.exit(1)

    outfile.write('{sample:{width}}'.format(sample='Sample', width=col_size(results)))

    header_elems = ['Date','MAPD','RNA_Reads','Expr_Sum']

    if dna_only:
        header_elems = header_elems[0:2]
        fstring = '{:<14}{:<10}\n'
    else:
        fstring = '{:<14}{:<10}{:<14}{:<14}\n' 
        if l[0] == 6:
            fstring = fstring.replace('\n','{:<14}{:<14}\n')
            header_elems += ['Pool1','Pool2']

    outfile.write(fstring.format(*header_elems))

    for sample in results:
        outfile.write('{sample:{width}}'.format(sample=sample,width=col_size(results)))
        out_res = [results[sample][r] for r in header_elems]
        outfile.write(fstring.format(*out_res))

def main():
    args = get_args()
    out_fh=''
    if args.output:
        sys.stdout.write('Writing results to %s.\n' % args.output)
        out_fh = open(args.output,'w')
    else:
        out_fh = sys.stdout

    results = {}
    for vcf in args.vcf:
        sample_name = get_name_from_vcf(vcf)
        results[sample_name] = read_vcf(vcf)
    print_data(results,out_fh,args.dna_only)

if __name__=='__main__':
    main()
