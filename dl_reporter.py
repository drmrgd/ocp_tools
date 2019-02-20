#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Read an Oncomine VCF and determine if there are outside labs variants to 
# report
# 
# 2019.02.11 - D Sims
################################################################################
"""
Starting with an Ion Reporter VCF that has been run as a part of the Oncomine
Comprehensive Assay (OCA) system, run through the MATCH MOI Reporter script and 
generate a list of variants that would be inclusionary for the NCI-MATCH Outside
Labs (AKA Designated Labs) initiative.
"""
import os
import sys
import csv
import argparse
import subprocess

from pprint import pprint as pp

from matchbox_api_utils import TreatmentArms

version = '0.6.022019'
match_arms = TreatmentArms(matchbox='adult', quiet=True)
sys.stderr.write("Note: using version %s of Treatment Arms DB.\n" 
    % match_arms.db_date)

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'vcf', 
        metavar="<VCF File>",
        help = 'Input VCF file to run.'
    )
    parser.add_argument(
        '-o', '--outfile', 
        #default='output.csv', 
        metavar='<outfile>',
        help='Custom output file (DEFAULT: %(default)s)'
    )
    parser.add_argument(
        '-s', '--status',
        metavar = "<Arm Status>",
        choices = ['ALL', 'OPEN', 'CLOSED', 'SUSPENDED'],
        default = 'OPEN',
        help='Only output arms with this status. Default: %(default)s.'
    )
    parser.add_argument(
        '-d', '--dl_excluded',
        action="store_false",
        help='Do not restrict output to arms that are open to the Designated '
            'Labs program.'
    )
    parser.add_argument(
        '-v', '--version', 
        action='version', 
        version = '%(prog)s - ' + version
    )
    args = parser.parse_args()

    ol_string = "are" if args.dl_excluded else "are not"
    sys.stderr.write("Outputting arms with status %s that %s open to outside "
        "labs.\n" % (args.status, ol_string))
    return args

def read_vcf(vcf, status, outside):
    """
    Run the VCF file through match_moi_report.pl and generate a CSV of data.
    """
    sys.stdout.write("Running `match_moi_report.pl` on VCF file...\n")
    cmd = ['match_moi_report.pl', '--cn' , '7', '--reads', '1000', '--Raw', vcf]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
        encoding="UTF-8")
    pout, perr = p.communicate()
    if p.returncode != 0:
        sys.stderr.write('ERROR: Encountered problem running '
            '`match_moi_report.pl`. Returned error:\n')
        sys.stderr.write("%s" % perr)
        sys.stderr.flush()
    else:
        var_data = csv.reader(pout.split('\n'), delimiter=',')

        # Always get a blank line in the output, which we don't want. 
        return build_variant_dict(list(var_data)[:-1], status, outside)

def build_variant_dict(variant_data, status, outside):
    """
    Read each variant line and build the required dictionary for parsing with 
    matchbox_api_utils.
    """
    for var in variant_data:
        # Set up a dict to pass into the amoi mapper function. Need all keys,
        # so set unecessary values to `None`.
        var_query = dict((x, None) for x in ('type', 'gene', 'identifier',
            'exon', 'function', 'oncominevariantclass'))
        if var[0] == 'SNV':
            var_query['type']                 = 'snvs_indels'
            var_query['gene']                 = var[9]
            var_query['identifier']           = var[8]
            var_query['exon']                 = var[13] 
            var_query['function']             = var[14]
            var_query['oncominevariantclass'] = var[15]
        elif var[0] == 'CNV':
            var_query['type'] = 'cnvs'
            var_query['gene'] = var[1]
        elif var[0] == 'Fusion':
            var_query['type']       = 'fusions'
            var_query['gene']       = var[4]
            var_query['identifier'] = '{}.{}'.format(var[1], var[2])

        if status == 'ALL':
            status = None

        # Add the aMOI mapping data from MATCHbox.
        arms = match_arms.map_amoi(var_query, outside=outside, status=status)
        if arms:
            var.append(';'.join(arms))
        else:
            var.append('---')
    return variant_data

def print_data(data, outfile):
    """
    Print the results out as a CSV with non-conforming title strings and whatnot
    to make somewhat useable in both CSV and human readable format.  If we start
    to use this programmatically, we can make it conventional CSV.
    """
    if outfile:
        sys.stderr.write('Writing data to %s\n' % outfile)
        outfh = open(outfile, 'w')
    else:
        outfh = sys.stdout
    csv_writer = csv.writer(outfh, delimiter=',', lineterminator='\n')

    snv_results = [var for var in data if var[0] == 'SNV']
    cnv_results = [var for var in data if var[0] == 'CNV']
    fusion_results = [var for var in data if var[0] == 'Fusion']
    
    outfh.write(':::  SNV Results :::\n')
    csv_writer.writerow(['Chr:Position', 'REF', 'ALT', 'VAF', 'Cov', 'ID',
        'Gene', 'Transcript', 'CDS', 'AA', 'Exon', 'Function', 'VariantClass',
        'MATCH_Arms'])
    if snv_results:
        fields = (1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 16, 17)
        for var in snv_results:
            csv_writer.writerow([var[i] for i in fields])
    else:
        outfh.write('No SNVs found.\n')
    outfh.write('\n')

    outfh.write(':::  CNV Results :::\n')
    csv_writer.writerow(['Chr', 'Gene', 'CN', 'MAPD', 'MATCH_Arms'])
    if cnv_results:
        fields = (2, 1, 5, 7, 8)
        for var in cnv_results:
            csv_writer.writerow([var[i] for i in fields])
    else:
        outfh.write('No CNVs found.\n')
    outfh.write('\n')

    outfh.write(':::  Fusion Results :::\n')
    csv_writer.writerow(['Fusion', 'ID', 'Drive_Gene', 'Reads', 'MATCH_Arms'])
    if fusion_results:
        fields = (1, 2, 4, 3, 6)
        for var in fusion_results:
            csv_writer.writerow([var[i] for i in fields])
    else:
        outfh.write('No Fusions found.\n')
    outfh.write('\n')

def main(vcf, outfile, arm_status, dl_excluded):
    variant_data = read_vcf(vcf, arm_status, dl_excluded)
    print_data(variant_data, outfile)

if __name__ == '__main__':
    args = get_args()
    main(args.vcf, args.outfile, args.status, args.dl_excluded)
