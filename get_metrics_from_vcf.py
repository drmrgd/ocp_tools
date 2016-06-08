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
from pprint import pprint


version = '2.0.0_060816'

def read_vcf(vcf_file):
    mapd_value = ''
    tot_rna_reads = ''
    expr_sum = 0

    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith('##mapd='):
                mapd_value = get_value(line.rstrip())
            elif line.startswith('##TotalMappedFusionPanelReads='):
                elems = line.split('=')
                tot_rna_reads = value = get_value(line.rstrip())
            elif re.search('SVTYPE=ExprControl',line):
                read_count = re.search('READ_COUNT=(\d+).*',line).group(1)
                expr_sum += int(read_count)
    return (mapd_value,tot_rna_reads,expr_sum)


def get_value(line):
    return line.split('=')[1]

def get_name(vcf):
    name_elems = vcf.split('_')
    sample_name = ''
    dna_name = ''
    rna_name = ''

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
            # sys.stderr.write("WARN: Can not determine DNA sample name from VCF filename. Using filename instead\n")
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

def print_data(results):
    print '{sample:{width}}{mapd:8}{rna_reads:12}{expr}'.format(
            sample='Sample', width=col_size(results), mapd='MAPD',rna_reads='RNA_Reads',expr='Expr_Sum')
    for sample in results:
        print '{sample:{width}}{mapd:<8}{rna_reads:<12}{expr}'.format(
                sample=sample,width=col_size(results),mapd=results[sample][0],rna_reads=results[sample][1],expr=results[sample][2])
    
def main():
    try:
        vcf_files = sys.argv[1:]
    except IndexError:
        sys.stderr.write("ERROR: No VCF files loaded to script!  You must enter at least one VCF file!\n")
        sys.exit(1)

    results = {}
    for vcf in vcf_files:
        sample_name = get_name(vcf)
        (results[sample_name]) = read_vcf(vcf)
    print_data(results)
    return

if __name__=='__main__':
    main()
