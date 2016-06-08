#!/usr/bin/python
# Since I have to frequently grab this table, output a table of MAPD values for a particular VCF or set of 
# VCF files.
#
# 3/28/2016 - D Sims
###########################################################################################################
import sys
import os
import re
from pprint import pprint

def read_vcf(vcf_file):
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith('##mapd='):
                elems = line.split('=')
                return elems[1].rstrip()

def get_name(vcf):
    sample_name = ''
    name_elems = vcf.split('_')
    try: 
        name_elems[0] == name_elems[2]
        sample_name = name_elems[0]
    except IndexError:
        try:
            dna_string = re.search(r'(.*?)-DNA_.*',vcf)
            sample_name = dna_string.group(1)
        except:
            sys.stderr.write("WARN: Can not determine DNA sample name from VCF filename. Using filename instead\n")
            sample_name = vcf.rstrip('.vcf')
    return sample_name

def col_size(data):
    col_width = 0
    for i in data:
        if len(i) > col_width:
            col_width = len(i)
    return col_width + 4

def print_data(results,delimiter):
    print '{sample:{width}}{delimiter}{mapd}'.format(
            sample='Sample', width=col_size(results), delimiter=delimiter, mapd='MAPD')
    for sample in results:
        print '{sample:{width}}{delimiter}{mapd}'.format(
                sample=sample,width=col_size(results),delimiter=delimiter,mapd=results[sample])
    
def main():
    try:
        vcf_files = sys.argv[1:]
    except IndexError:
        sys.stderr.write("ERROR: No VCF files loaded to script!  You must enter at least one VCF file!\n")
        sys.exit(1)

    results = {}
    for vcf in vcf_files:
        sample_name = get_name(vcf)
        results[sample_name] = read_vcf(vcf)

    print_data(results,' ')
    return

if __name__=='__main__':
    main()
