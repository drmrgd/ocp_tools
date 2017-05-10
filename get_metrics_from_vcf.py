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

version = '3.3.2_031017'

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
                    fetched_data['MAPD'] = (get_value(line.rstrip()))

                elif line.startswith('##TotalMappedFusionPanelReads='):
                    fetched_data['RNA_Reads'] = (get_value(line.rstrip()))

                elif line.startswith('##fileDate'):
                    date = datetime.datetime.strptime(get_value(line.rstrip()),"%Y%M%d")
                    formatted_date = date.strftime("%Y-%M-%d")
                    fetched_data['Date'] = (formatted_date)

                elif line.startswith('##OncomineVariantAnnotationToolVersion'):
                    ovat_version = get_value(line.rstrip())
                    if LooseVersion(ovat_version) > LooseVersion(oca_v3_version):
                        p1,p2 = get_rna_pool_info(vcf_file)
                        fetched_data['Pool1'] = p1
                        fetched_data['Pool2'] = p2

                elif re.search('SVTYPE=ExprControl',line):
                    read_count = re.search('READ_COUNT=(\d+).*',line).group(1)
                    expr_sum += int(read_count)
            fetched_data['Expr_Sum'] = str(expr_sum)
    except IOError as e:
        sys.stderr.write('ERROR: Can not open file {}: {}!\n'.format(vcf_file,e))
        sys.exit()
    return fetched_data

def get_rna_pool_info(vcf):
    p = subprocess.Popen(['match_rna_qc.pl', vcf], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (data,err) = p.communicate()
    ret_res = data.split('\n')
    results = dict(zip(ret_res[0].split(','),ret_res[1].split(',')))
    # return int(round(float(results['pool1_total']))), int(round(float(results['pool2_total'])))
    p1_tot = str(round(float(results['pool1_total'])))
    p2_tot = str(round(float(results['pool2_total'])))
    return p1_tot.split('.')[0], p2_tot.split('.')[0]

def get_value(line):
    return line.split('=')[1]

def get_name(vcf):
    name_elems = vcf.split('_')
    sample_name = ''
    dna_name = ''
    rna_name = ''

    # If this is MATCHBox data, we always start with 'MSN####'
    if name_elems[0].startswith('MSN'):
        return name_elems[0]

    try:
        if name_elems[0] == name_elems[2]:
            dna_name = name_elems[0]
            rna_name = name_elems[2]
    except IndexError:
        try:
            (dna_name, rna_name) = re.search(r'^(.*?(?:DNA|_v\d)?)_(.*?RNA).*$',vcf).group(1,2)
            if not dna_name.endswith('DNA'):
                dna_name += '-DNA'
        except:
            #sys.stderr.write("WARN: Can not determine DNA sample name from VCF filename. Using filename instead\n")
            vcf_name = vcf.rstrip('.vcf')

    if dna_name and rna_name:
        sample_name = '_'.join([dna_name,rna_name])
    else:
        sample_name = vcf_name
    return (sample_name) 

def col_size(data):
    col_width = 0
    for i in data:
        if len(i) > col_width:
            col_width = len(i)
    return col_width + 4

def print_data(results,outfile):
    # Figure out if we have two different versions of analysis, and if so bail out to make easier.
    l = [len(v) for k,v in results.iteritems()]
    if len(set(l)) > 1:
        print('Mixed version VCFs detected!  We can not process two different versions together!  Please run separately and cat the data later.')
        sys.exit(1)

    header_elems = ['Date','MAPD','RNA_Reads','Expr_Sum']
    fstring = '{:12}{:8}{:12}{:12}\n' 
    outfile.write('{sample:{width}}'.format(sample='Sample', width=col_size(results)))

    if l[0] == 6:
        fstring = fstring.replace('\n','{:12}{:12}\n')
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
        sample_name = get_name(vcf)
        results[sample_name] = read_vcf(vcf)

    print_data(results,out_fh)

if __name__=='__main__':
    main()
