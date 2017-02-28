#!/usr/bin/python
# Run MATCH MOI Reporter on a VCF file and generate a CSV with result that can be imported into Excel for a collated
# report
#
# 2/11/2014
#######################################################################################################################
import sys
import os
import re
import subprocess
import argparse
from natsort import natsorted
from collections import defaultdict
from pprint import pprint
from multiprocessing.pool import ThreadPool

version = '2.1.1_021517'

def get_args():
    parser = argparse.ArgumentParser(
            formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 100, width=200),
            description = '''
            Collect MOI reports for a set of VCF files and output a raw CSV formatted output that can be
            easily imported into Microsoft Excel for reporting.  This script relies on `match_moi_report.pl`
            in order to generate the MOI reports for each.
            ''',
            version = '%(prog)s - ' + version,
            )
    parser.add_argument('vcf_files', nargs='+', help='List of VCF files to process.')
    parser.add_argument('-q','--quiet', action='store_false', default=True, 
            help='Do not suppress warning and extra output')
    parser.add_argument('--cu', default=4, metavar='INT', 
            help='Copy number threshold (5%% CI lower bound) to pass to match_moi_report for reporting amplifications. DEFAULT: %(default)s')
    parser.add_argument('--cl', default=1, metavar='INT', 
            help='Copy number threshold (95%% CI upper bound) to pass to match_moi_report for reporting copy loss. DEFAULT: %(default)s')
    parser.add_argument('--reads', default=100, metavar='INT',
            help='Threshold for number of fusion reads to report. DEFAULT: %(default)s')
    parser.add_argument('-o','--output', help='Output to file rather than STDOUT ***NOT YET IMPLEMENTED***')
    args = parser.parse_args()

    global quiet
    quiet = args.quiet
    return args

def get_names(string):
    string = os.path.basename(string)
    match_1 = re.search('(.*?)[-_]?DNA_(.*?)[-_]?RNA.*\.vcf', string)    
    try: 
        dna_samp = match.group(1)
        rna_samp = match.group(2)
    except:
        if not quiet:
            sys.stderr.write("WARN: Can not get DNA or RNA sample name for '%s'! Using full VCF filename instead\n" % string)
            
        dna_samp = rna_samp = string.rstrip('.vcf')
    return dna_samp, rna_samp

# def gen_moi_report(vcf,dna,rna):
def gen_moi_report(vcf,cu,cl,reads):
    '''Use MATCH MOI Reporter to generate a variant table we can parse later. Gen CLI Opts to determine
    what params to run match_moi_report with'''
    (dna,rna) = get_names(vcf)
    moi_report_cmd = ['match_moi_report.pl', '--cu', str(cu), '--cl', str(cl), '-r', str(reads), '-R', vcf]
    p=subprocess.Popen(moi_report_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result,error = p.communicate()

    if p.returncode != 0:
        sys.stderr.write("ERROR: Can not process file: {}!\n".format(vcf))
        raise(error)
    else:
        return vcf, parse_data(result,dna,rna)

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

def print_data(var_type,data,outfile):
    # Split the key by a colon and sort based on chr and then pos using the natsort library
    if var_type == 'null':
        # print ','.join(data['no_result'])
        outfile.write(','.join(data['no_result']) + "\n")
    elif var_type == 'snv_data':
        for variant in natsorted(data.keys(), key=lambda k: (k.split(':')[1], k.split(':')[2])):
            # print ','.join(data[variant])
            outfile.write(','.join(data[variant]) + "\n")
    else:
        for variant in natsorted(data.keys(), key=lambda k: k.split(':')[1]):
            # print ','.join(data[variant])
            outfile.write(','.join(data[variant]) + "\n")
    return

def arg_star(args):
    return gen_moi_report(*args)

def main():
    args = get_args()

    # Single process
    # moi_data = defaultdict(dict)
    # for x in args.vcf_files:
        # print 'processing %s...' %  x
        # moi_data[x] = gen_moi_report(x,args.cu,args.cl,args.reads)
    # pprint(dict(moi_data))
    # sys.exit()

    # Threaded method
    pool = ThreadPool(48)
    task_list = [(x,args.cu,args.cl,args.reads) for x in args.vcf_files]

    try:
        moi_data = {vcf : data for vcf,data in pool.imap_unordered(arg_star,task_list)}
    except Exception:
        pool.close()
        pool.join()
        sys.exit(1)

    # Setup an output file if we want one
    outfile = ''
    if args.output:
        print "Writing output to '%s'" % args.output
        outfile = open(args.output, 'w')
    else:
        outfile = sys.stdout

    header = ['Sample','Gene','Position','Ref','Alt','VARID','Type','VAF/CN','Coverage/Counts','RefCov','AltCov']
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
