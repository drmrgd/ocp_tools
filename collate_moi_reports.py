#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Run MATCH MOI Reporter on a VCF file and generate a CSV with result that can 
# be imported into Excel for a collated report
#
# 2/11/2014
################################################################################
"""
Collect MOI reports for a set of VCF files and output a raw CSV formatted output 
that can be easily imported into Microsoft Excel for reporting. This script relies 
on `match_moi_report.pl` in order to generate the MOI reports for each.
"""
import sys
import os
import re
import subprocess
import argparse
from natsort import natsorted
from collections import defaultdict
from pprint import pprint as pp
from multiprocessing.pool import ThreadPool

version = '2.11.012918'
debug = False 

def get_args():
    # Default thresholds. Put them here rather than fishing below.
    num_procs = '24'
    cn = 7
    cu = None
    cl = None
    reads = 1000

    parser = argparse.ArgumentParser(
            formatter_class = lambda prog: argparse.HelpFormatter(prog, 
                max_help_position = 100, width=200),
            description = __doc__,
    )
    parser.add_argument('vcf_files', nargs='+', help='List of VCF files to process.')
    parser.add_argument('-j', '--procs', default=num_procs, metavar='<num_procs>',
        help='Number of processes to run. Input 0 if you want do not want to perform '
            'threaded processing. DEFAULT: %(default)s')
    parser.add_argument('-q','--quiet', action='store_false', default=True, 
        help='Do not suppress warning and extra output')
    parser.add_argument('--cn', metavar='INT', default = cn,
        help='Use copy number (CN) value for CNV reporting to be more compatible '
            'with MATCH rules. This will disable CU and CL thresholds and is on '
            'by default')
    parser.add_argument('--cu', default=cu, metavar='INT', 
        help='Copy number threshold (5%% CI lower bound) to pass to match_moi_report '
            'for reporting amplifications. DEFAULT: %(default)s')
    parser.add_argument('--cl', default=cl, metavar='INT',
        help='Copy number threshold (95%% CI upper bound) to pass to match_moi_report '
            'for reporting copy loss. DEFAULT: %(default)s')
    parser.add_argument('--reads', default=reads, metavar='INT', 
        help='Threshold for number of fusion reads to report. DEFAULT: %(default)s')
    parser.add_argument('-o','--output', help='Output to file rather than STDOUT.')
    parser.add_argument('-p', '--pedmatch', action='store_true', 
        help='Data comes from Pediatric MATCH rather than Adult MATCH (default).')
    parser.add_argument('-b', '--blood', action='store_true', 
        help='Data comes from blood specimens, and therefore we only have DNA data.')
    parser.add_argument('-v', '--version', action='version', 
        version = '%(prog)s - ' + version)
    args = parser.parse_args()

    if args.cn:
        args.cu = args.cl = None

    global quiet
    quiet = args.quiet
    return args

def get_names(string):
    string = os.path.basename(string)
    try:
        (dna_samp,rna_samp) = re.search('(.*?DNA)_(.*).vcf$',string).group(1,2)
    except:
        if not quiet:
            sys.stderr.write("WARN: Can not get DNA or RNA sample name for '%s'! "
                "Using full VCF filename instead\n" % string)
        dna_samp = rna_samp = string.rstrip('.vcf')
    return dna_samp, rna_samp

def parse_cnv_params(cu,cl,cn):
    '''
    Since CNV params are a bit difficult to work with, create a better, standardized 
    list to pass into functions below
    '''
    params_list = []
    params = {
        '--cu' : cu,
        '--cl' : cl,
        '--cn' : cn 
    }

    # Would like to write a generator expression here, but can't figure it out!
    for k,v in params.iteritems():
        if v:
            params_list.extend([k,str(v)])
    return params_list

def gen_moi_report(vcf,cnv_args,reads,proc_type,study,blood):
    '''
    Use MATCH MOI Reporter to generate a variant table we can parse later. Gen 
    CLI Opts to determine what params to run match_moi_report with
    '''
    (dna,rna) = get_names(vcf)

    # thresholds = cnv_args + ['-r', str(reads), '-R']
    thresholds = cnv_args

    # Determine if blood or tumor and add appropriate thresholds
    if not blood:
        thresholds += ['-r', str(reads), '-R']
    else:
        thresholds.append('-Rb')

    # If Pediatric MATCH need to pass different MOI rules.
    if study:
        thresholds.append('-p')

    moi_report_cmd = ['match_moi_report.pl'] + thresholds + [vcf]
    p=subprocess.Popen(moi_report_cmd, stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE)
    result,error = p.communicate()

    if p.returncode != 0:
        sys.stderr.write("ERROR: Can not process file: {}!\n".format(vcf))
        raise(error)
    else:
        # need a tuple to track threads and not crash dict entries if we're 
        # doing multithreaded processing.
        if proc_type == 'single':
            return parse_data(result,dna,rna)
        elif proc_type == 'threaded':
            return vcf, parse_data(result,dna,rna)

def populate_list(var_type, var_data):
    wanted_fields = {
        'snv'     : [0,9,1,2,3,10,11,12,8,4,5,6,7,13],
        'cnv'     : [0,1,2,5],
        'fusions' : [0,4,2,1,3]
    }
    return [var_data[x] for x in wanted_fields[var_type]]

def pad_list(data_list, data_type):
    '''
    Pad out the list with hyphens where there is no relevent data.  Maybe kludgy, 
    but I don't know a better way
    '''
    tmp_list = ['-'] * 13
    data_list.reverse()
    if data_type == 'cnv':
        for i in [0,1,2,9]:
            tmp_list[i] = data_list.pop()
    elif data_type == 'fusions':
        # Get rid of 
        for i in [0,1,8,4,10]:
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
        data['null']['no_result'] = [dna] + ['-']*11
    return dict(data)

def print_data(var_type,data,outfile):
    # Split the key by a colon and sort based on chr and then pos using the 
    # natsort library
    if var_type == 'null':
        outfile.write(','.join(data['no_result']) + "\n")
    elif var_type == 'snv_data':
        for variant in natsorted(
                data.keys(), key=lambda k: (k.split(':')[1], k.split(':')[2])):
            outfile.write(','.join(data[variant]) + "\n")
    else:
        for variant in natsorted(data.keys(), key=lambda k: k.split(':')[1]):
            outfile.write(','.join(data[variant]) + "\n")
    return

def arg_star(args):
    return gen_moi_report(*args)

def non_threaded_proc(vcf_files,cnv_args,reads,study,blood):
    '''
    If for some reason we only want to use a single thread / proc to do this
    '''
    moi_data = defaultdict(dict)
    for x in vcf_files:
        print 'processing %s...' %  x
        moi_data[x] = gen_moi_report(x,cnv_args,reads,'single',study,blood)
    return moi_data

def threaded_proc(vcf_files,cnv_params,reads,study,blood):
    pool = ThreadPool(48)
    moi_data = defaultdict(dict)
    task_list = [(x,cnv_params,str(reads),'threaded',study,blood) for x in vcf_files]
    try:
        moi_data = {vcf : data for vcf, data in pool.imap_unordered(arg_star, task_list)}
    except Exception:
        raise
        pool.close()
        pool.join()
        sys.exit(1)
    return moi_data

def print_title(fh,**kwargs):
    '''Print out a header to remind me just what params I used this time!'''
    cnv_params = parse_cnv_params(kwargs['cu'],kwargs['cl'],kwargs['cn'])
    string_params = '='.join([str(x).lstrip('--') for x in cnv_params])
    if 'cl' in string_params:
        string_params = string_params.replace('=cl','; cl')
    if kwargs['pedmatch']:
        study_name = 'Pediatric MATCH'
    else:
        study_name = 'Adult MATCH'

    fh.write('-'*95)
    fh.write('\nCollated {} MOI Reports Using Params CNV: {}, Fusion Reads: '
        'reads={}\n'.format(study_name,string_params,kwargs['reads']))
    fh.write('-'*95)
    fh.write('\n')
    return

def main():
    args = get_args()
    if debug:
        print('CLI args as passed:')
        pp(vars(args))
        print('')

    # Setup an output file if we want one
    outfile = ''
    if args.output:
        print "Writing output to '%s'" % args.output
        outfile = open(args.output, 'w')
    else:
        outfile = sys.stdout

    # handle complex CNV args
    cnv_args = parse_cnv_params(args.cu,args.cl,args.cn)

    if args.procs == '0':
        moi_data = non_threaded_proc(args.vcf_files, cnv_args, args.reads,
            args.pedmatch, args.blood)
    else:
        moi_data = threaded_proc(args.vcf_files, cnv_args, args.reads, 
            args.pedmatch, args.blood)

    # Print data
    print_title(outfile,**vars(args))
    header = ['Sample', 'Type', 'Gene', 'Position', 'Ref', 'Alt', 
        'Transcript', 'CDS', 'AA', 'VARID', 'VAF/CN', 'Coverage/Counts',
        'RefCov', 'AltCov', 'Function']
    outfile.write(','.join(header) + "\n")
    
    # Print out sample data by VCF
    var_types = ['snv_data', 'cnv_data', 'fusion_data', 'null']
    for sample in sorted(moi_data):
        for var_type in var_types:
            try:
                print_data(var_type,moi_data[sample][var_type],outfile)
            except KeyError:
                continue

if __name__ == '__main__':
    main()
