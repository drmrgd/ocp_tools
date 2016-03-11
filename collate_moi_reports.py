#!/usr/bin/python
# Run MATCH MOI Reporter on a VCF file and generate a CSV with result that can be imported into Excel for a collated
# report
#
# TODO:
#    - Add some CLI opts for output options and args to pass to match_moi_report.pl
#    - Can we parallelize the runs to speed them up through the MOI reporter?
#
# 2/11/2014
#######################################################################################################################
import sys
import os
import re
import subprocess
import operator
from natsort import natsorted
from collections import defaultdict
from pprint import pprint

VERSION = '0.9.2_022316'

def get_names(string):
    string = os.path.basename(string)
    match = re.search('(.*?)[-_]?DNA_(.*?)[-_]?RNA.*\.vcf', string)    
    try: 
        dna_samp = match.group(1)
        rna_samp = match.group(2)
    except:
        sys.stderr.write("WARN: Can not get DNA or RNA sample name for '%s'! Using full VCF filename instead\n" % string)
        dna_samp = rna_samp = string.rstrip('.vcf')
    return dna_samp, rna_samp

def gen_moi_report(vcf,dna,rna):
    '''Use MATCH MOI Reporter to generate a variant table we can parse later. Gen CLI Opts to determine
    what params to run match_moi_report with'''
    p=subprocess.Popen(['match_moi_report.pl','-r','-c4','-f0.03',vcf], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #p=subprocess.Popen(['match_moi_report.pl','-r',vcf], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result,error = p.communicate()
    if p.returncode != 0:
        sys.stderr.write("ERROR: Can not process file: {}!".format(vcf))
        sys.exit(1)
    else:
        return parse_data(result,dna,rna)

def populate_list(var_type, var_data):
    wanted_fields = {
        'snv'     : [9,1,2,3,8,0,4,5,6,7],
        'cnv'     : [1,2,0,5],
        'fusions' : [4,2,1,0,3]
    }
    return [var_data[x] for x in wanted_fields[var_type]]

def pad_list(data_list,data_type):
    '''Pad out the list with hyphens where there is no relevent data.  Maybe kludgy, but I don't know a better way'''
    tmp_list = ['-'] * 10
    data_list.reverse()
    if data_type == 'cnv':
        for i in [0,1,5,6]:
            tmp_list[i] = data_list.pop()
    elif data_type == 'fusions':
        for i in [0,3,4,5,7]:
            tmp_list[i] = data_list.pop()
    return tmp_list

def parse_data(report_data,dna,rna):
    data = defaultdict(dict)
    raw_data = report_data.split('\n')

    for line in raw_data:
        fields = line.split(',')
        if fields[0] == 'SNV':
            varid = fields[9] +':'+ fields[1]
            data['snv_data'][varid] = [dna] + populate_list('snv', fields)
        elif fields[0] == 'CNV':
            varid = fields[1] +':'+ fields[2]
            cnv_data = populate_list('cnv', fields)
            padded_list = pad_list(cnv_data,'cnv')
            data['cnv_data'][varid] = [dna] + padded_list
        elif fields[0] == 'Fusion':
            varid = fields[1] +':'+ fields[2]
            fusion_data = populate_list('fusions', fields)
            padded_list = pad_list(fusion_data,'fusions')
            data['fusion_data'][varid] = [rna] + padded_list

    # Let's still output something even if no MOIs were detected
    if not data:
        data['null']['no_result'] = [dna] + ['-']*10
    return data

def print_data(var_type,data):
    # Split the key by a colon and sort based on chr and then pos using the natsort library
    if var_type == 'snv_data':
        for variant in natsorted(data.keys(), key=lambda k: (k.split(':')[1], k.split(':')[2])):
            # sys.stdout.write('{} => '.format(variant))
            print ','.join(data[variant])
    else:
        for variant in natsorted(data.keys(), key=lambda k: k.split(':')[1]):
            print ','.join(data[variant])
    return

def main():
    vcf_files = sys.argv[1:]
    moi_data = defaultdict(dict)

    for vcf in vcf_files:
        (dna, rna) = get_names(vcf)
        # print "dna: {}\nrna: {}".format(dna, rna)
        moi_data[vcf] = gen_moi_report(vcf,dna,rna)

    # Set up and print report header
    header = ['Sample', 'Gene', 'Position', 'Ref', 'Alt', 'VARID', 'Type', 'VAF/CN', 'Coverage/Counts', 
            'RefCov', 'AltCov']
    print ','.join(header)
    
    # Print out sample data by VCF
    var_types = ['snv_data', 'cnv_data', 'fusion_data']
    for sample in sorted(moi_data):
        for var_type in var_types:
            print_data(var_type,moi_data[sample][var_type])

if __name__ == '__main__':
    main()
