#!/usr/bin/python
import sys
import os
import re
import subprocess
import argparse
from natsort import natsorted
from pprint import pprint
from collections import defaultdict

version = '1.3.0_022916'

def get_opts():
    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 100, width=200),
        description = '''
        Input a list of NCI-MATCH VCF files from MATCHBox and output a collated CNV dataset for a specific gene.  Must 
        have 'ocp_cnv_report.pl' in your path in order to work.
        ''',
        version = '%(prog)s - ' + version,
        )
    parser.add_argument('gene', metavar='gene_name', help='Gene name to query. NOTE: name is case sensitive')
    parser.add_argument('vcf_files', nargs='+', metavar='vcf_files', help='VCF files to process')
    parser.add_argument('-n', '--nonmatch', action='store_true', 
            help='VCF files are not from MATCHBox and do not contain PSN / MSN information')
    parser.add_argument('-cn', '--copy-number', metavar='INT', help='Only print results if CN is greater that this value')
    parser.add_argument('-csv', action='store_true', 
            help='Format results as a CSV file (default is whitespace formatted, pretty print file')
    parser.add_argument('-tsv', action='store_true', 
            help='Format results as a TSV file (default is whitespace formatted, pretty print file')
    args = parser.parse_args()

    return args

def read_file(vcf,gene,cn):
    cmd = 'ocp_cnv_report.pl -g {} -c {} {}'.format(gene, str(cn), vcf)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result,err = p.communicate()
    try:
        cnv_data = parse_raw_data(result)
    except:
        print "WARNING: error processing file {}. Skipping this VCF file.".format(vcf)
        return
    return cnv_data

def convert_type(i):
    data_types = [int, float, str]
    for t in data_types:
        try:
            return t(i)
        except ValueError:
            pass

def parse_raw_data(raw_data):
    lines = raw_data.split('\n')
    data = []

    for line in lines:
        if line.startswith('chr'):
            elems = line.split()
            data.extend(line.split())
        elif line.startswith(':::'):
            match = re.match('.*Data For (\w+).*Gender: (Unknown|Male|Female).*MAPD: (.*)\) :::$',line)
            sample_name = match.group(1)
            gender = match.group(2)
            mapd = match.group(3)
            data.extend([sample_name,gender,mapd])

    # pad out with "NDs" if no match.
    if len(data) < 4:
        null_result = 'ND '*12
        data.extend([x for x in null_result.split()])

    clist = [convert_type(x) for x in data]
    for i in range(9,14,1):
        if type(clist[i]) != str:
            clist[i] = format(clist[i], '.1f')
    return clist

def output_data(cnv_data,delimiter,nonmatch):
    template = ''

    if delimiter:
        template = ('{}'+delimiter)*8 + '{}'
    else:
        template = '{:<11} {:<10} {:<10} {:<8} {:<10} {:<8} {:<7} {:<7} {:<3}'

    if nonmatch:
        header = ('Patient', 'MSN', 'Gender', 'MAPD', 'Gene', 'Chr', 'CI_05', 'CI_95', 'CN')
    else:
        header = ('Filename', 'Sample', 'Gender', 'MAPD', 'Gene', 'Chr', 'CI_05', 'CI_95', 'CN')
    print template.format(*header)

    desired_elements = [0,1,2,4,3,11,12,13]
    for patient in natsorted(cnv_data):
        print template.format(patient,*[cnv_data[patient][i] for i in desired_elements])
    return

def main():
    args = get_opts()
    gene_lookup = args.gene
    cn_threshold = 0
    if args.copy_number:
        cn_threshold = args.copy_number

    delimiter = None
    if args.csv:
        delimiter = ','
    elif args.tsv:
        delimiter = '\t'

    results = defaultdict(list)
    for vcf in args.vcf_files:
        if args.nonmatch:
            (psn,msn,ver) = os.path.basename(vcf).split('_')
            ver=ver.rstrip('.vcf')
            results[psn] = read_file(vcf, gene_lookup, cn_threshold)
        else:
            results[vcf.rstrip('.vcf')] = read_file(vcf, gene_lookup, cn_threshold)

    output_data(results, delimiter, args.nonmatch)
    return

if __name__ == '__main__':
    main()
