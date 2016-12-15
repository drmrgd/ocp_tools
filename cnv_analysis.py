#!/usr/bin/python
import sys
import os
import re
import subprocess
import argparse
from natsort import natsorted
from pprint import pprint as pp
from collections import defaultdict
from multiprocessing.pool import ThreadPool

version = '2.0.1_000716'

def get_opts():
    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 100, width=200),
        description = '''
        Input a list of NCI-MATCH VCF files from MATCHBox and output a collated CNV dataset for a specific gene.  Must 
        have 'ocp_cnv_report.pl' in your path in order to work.
        ''',
        version = '%(prog)s - ' + version,
        )
    parser.add_argument('vcf_files', nargs='+', metavar='vcf_files', help='VCF files to process')
    parser.add_argument('-g', '--gene', metavar='gene_name', required=True,
            help='Gene name to query. Can input a comma separated list of genes to query.')
    parser.add_argument('-cn', '--copy-number', default=0, metavar='INT', help='Only print results if CN is greater that this value')
    parser.add_argument('-csv', action='store_true', 
            help='Format results as a CSV file (default is whitespace formatted, pretty print file')
    args = parser.parse_args()

    return args

def read_file(vcf,gene,cn):
    cmd = 'ocp_cnv_report.pl -g {} --cn {} {}'.format(gene, str(cn), vcf)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result,err = p.communicate()
    try:
        sample_name,cnv_data = parse_raw_data(result)
    except:
        print "WARNING: error processing file {}. Skipping this VCF file.".format(vcf)
        return
    return sample_name,cnv_data

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
    result = []

    for line in lines:
        if line.startswith('chr'):
            data.append(line.split())
        elif line.startswith(':::'):
            match = re.match('.*Data For (\w+).*Gender: (Unknown|Male|Female).*MAPD: (.*)\) :::$',line)
            sample_name = match.group(1)
            gender = match.group(2)
            mapd = match.group(3)
            sample_id = ':'.join([sample_name,gender,mapd])

    # pad out with "NDs" if no match.
    for elem in data: 
        if len(elem) < 4:
            null_result = 'ND '*12
            elem.extend([x for x in null_result.split()])

        clist = [convert_type(x) for x in elem]
        for i in range(6,11,1):
            if type(clist[i]) != str:
                clist[i] = format(clist[i], '.2f')
        result.append(clist)
    return sample_id,result

def output_data(cnv_data,delimiter):
    template = ''
    sample_width = get_width(cnv_data.keys())
    if delimiter:
        template = ('{}'+delimiter)*7 + '{}'
    else:
        template = '{:<{width}} {:<10} {:<8} {:<10} {:<8} {:<7} {:<7} {:<3}'
    header = ('Sample', 'Gender', 'MAPD', 'Gene', 'Chr', 'CI_05', 'CI_95', 'CN')
    print template.format(width=sample_width,*header)

    desired_elements = [1,0,8,9,10]
    for patient in natsorted(cnv_data):
        sample,gender,mapd = patient.split(':')
        for result in cnv_data[patient]:
            output_data = [sample,gender,mapd]
            output_data.extend([result[i] for i in desired_elements])
            print template.format(width=sample_width, *output_data)
    return

def arg_star(args):
    '''Constructor to generate a method call to read_file with supplied args.  We need this to be able to use the 
    pool function with multiple arguments. '''
    return read_file(*args)

def get_width(data):
    width = 0
    samples = [x.split(':')[0] for x in data]
    for sample in samples:
        if len(sample) > width:
            width = len(sample)
    return int(width) + 4

def main():
    args = get_opts()
    gene_lookup = args.gene
        
    cn_threshold = 0
    if args.copy_number:
        cn_threshold = args.copy_number

    delimiter = None
    if args.csv:
        delimiter = ','

    pool = ThreadPool(48)
    task_list = [(x,gene_lookup,cn_threshold) for x in args.vcf_files]
    results = {sample : data for sample,data in pool.imap_unordered(arg_star,task_list)}
    output_data(results, delimiter)
    return

if __name__ == '__main__':
    main()
