#!/usr/bin/python
import sys
import os
import re
import shutil
import argparse
import subprocess
import fnmatch
from time import sleep
from pprint import pprint

version = '1.6.1_121416'

def get_args():
    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 100, width=200),
        description='''
        Wrapper program for several shell commands and other programs that are invoked for a variant review process.
        Wrote a script to help automate this and (possibly) save on some carpal tunnel!
        ''',
        version = '%(prog)s  - ' + version,
        )
    parser.add_argument("dna_bam", help='DNA BAM file from MATCHBox.')
    parser.add_argument("rna_bam", help='RNA BAM file from MATCHBox.')
    parser.add_argument('-a', '--analysis_id', metavar='<ir_analysis_id>', 
            help='The IR analysis ID that will be used to pull the data from the IR server if different than in BAM name')
    parser.add_argument('-t', '--token', metavar='<ir_token_path>', 
            help='The path to the Ion Reporter token file used to pull the analysis from the server')
    parser.add_argument('-s', '--site', metavar='<match_site_id>', default='nci',
            help='The site ID used in the ir_api_retrieve script to pull out the data (DEFAULT: %(default)s)')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-p', '--psn', metavar='<PSN>', required=True,
            help='PSN for this case')
    requiredNamed.add_argument('-m', '--msn', metavar='<MSN>', required=True,
            help='MSN for this case')
    args = parser.parse_args()

    site = args.site
    psn  = args.psn
    msn  = args.msn 

    dna_bam = args.dna_bam
    if not os.path.isfile(dna_bam):
        sys.stderr.write("ERROR: '%s' does not exist!\n" % dna_bam)
        sys.exit(1)
        
    rna_bam = args.rna_bam
    if not os.path.isfile(rna_bam):
        sys.stderr.write("ERROR: '%s' does not exist!\n" % rna_bam)
        sys.exit(1)

    match=re.search('(PSN)([0-9]+$)', psn)
    if not match:
        sys.stderr.write("ERROR: '%s' is not a valid PSN!\n" % psn)
        sys.exit(1)

    match=re.search('(MSN)([0-9]+)$', msn)
    if not match:
        sys.stderr.write("ERROR: '%s' is not a valid MSN!\n" % msn)
        sys.exit(1)

    if args.analysis_id:
        run_id = args.analysis_id
    else:
        run_id = ''

    # Check for and hand the IR token file.
    token = os.environ['HOME'] + '/Dropbox/ir_stuff/ir_tokens'
    if args.token:
        token = args.token
    if not os.path.exists(token):
        sys.stderr.write("ERROR: The IR token file '%s' can not be found! You might need to load a custom token with the '-t' option.\n" % token)
        sys.exit(1)

    return (dna_bam, rna_bam, psn, msn, site, run_id, token)

def gen_wd(psn,msn):
    new_dir = psn +'_'+ msn +'_variant_reports'

    if os.path.isdir(new_dir):
        sys.stderr.write("WARN: Directory '%s' already exists!" % new_dir)
        choice = user_query(' Overwrite current data?')
        if not choice:
            sys.stdout.write("Exiting so we don't overwrite old data! You should move old data to a new directory and try again.\n")
            sys.exit(1)
        else:
            sys.stdout.write("Revoming old directory to make room for new data.\n")
            shutil.rmtree(new_dir)

    sys.stdout.write("Making new directory to store data: %s.\n" % new_dir)
    os.mkdir(new_dir)
    return new_dir

def user_query(query, default='no'):
    '''User query code from SO #3041986 user fmark.  Very nice idea!'''
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}

    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(query + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stderr.write("Invalid choice '%s'! Please enter 'y' or 'n'." % choice)

def validate_bams(msn, bam, na_type):
    '''Validate the BAM files passed into the script, and if OK, get a new file name, get analysis ID, and index them'''
    # match = re.search('^.*?(MSN\d+_.*?_[cv]\d+_.*)_([dr]na).bam', bam)
    # match = re.search('^.*?(MSN\d+_[cv]\d+_.*)_([dr]na).bam', bam)
    match = re.search('^.*?(MSN\d+_(?:[DR]NA_)?[cv]\d+_.*)_([dr]na).bam', bam)
    try:
        sample = match.group(2)
    except AttributeError:
        sys.stderr.write("ERROR: BAM file '%s' does not have valid MSN nomenclature and can not be processed.\n" % bam)
        sys.exit(1) 

    if sample != na_type:
        sys.stderr.write("ERROR: Expecting a %s file, but got file %s instead.  Check the file name and / or order!\n" % (na_type, bam))
        sys.exit(1)

    analysis_id = match.group(1)
    elems = analysis_id.split('_')
    if elems[0] != msn:
        sys.stderr.write("ERROR: the MSN '%s' does not agree with the MSN in the %s BAM file '%s'!\n" % (msn, na_type, bam))
        sys.exit(1)

    new_file_name = analysis_id + '_' + sample + '.bam'
    
    sys.stdout.write("\tStripping off extra prefix from %s file if exists..." % na_type.upper())
    os.rename(bam, new_file_name)
    sys.stdout.write("Done!\n")

    sys.stdout.write("\tIndexing %s BAM file..." % na_type.upper())
    p = subprocess.Popen(['samtools','index', new_file_name], stdout=subprocess.PIPE)
    p.communicate()
    sys.stdout.write("Done!\n")

    return (new_file_name, analysis_id)

def get_ir_data(analysis_id, token):
    '''Use ir_api_retrieve and extract_ir_data to get IR data for review'''

    sys.stdout.write('Getting data from IR for analysis ID {}...'.format(analysis_id))
    count = 1
    while count:
        p = subprocess.Popen(['ir_api_retrieve.py', '-H','nci', '-t', token, analysis_id],stdout=subprocess.PIPE,stderr=subprocess.PIPE) 
        result,error = p.communicate()
        if p.returncode != 0:
            sys.stderr.write("\nError retrieving data: {}\n".format(error))
            if count == 3:
                sys.stderr.write('\nERROR: Can not retrieve data from IR server. Max attempts reached (3).  Check the connection and try again!\n')
                sys.exit(1)
            sys.stderr.write('\nWaiting 5 seconds and trying again [attempt: {}]...\n'.format(count+1))
            sleep(5)
            count += 1
        else:
            break
    sys.stdout.write('Done!\n')

    sys.stdout.write("Extracting IR results...")
    p = subprocess.Popen('extract_ir_data.sh', stdout=subprocess.PIPE)
    p.communicate()
    sys.stdout.write("Done!\n")

def gen_moi_report(msn,vcf):
    filename = msn + '_MATCH_MOI_Report.txt'

    sys.stdout.write("Generating a MATCH MOI Report for {}...".format(msn))
    p=subprocess.Popen(['match_moi_report.pl', '--cu', '4', '--cl', '1', '-o', filename, vcf], stdout=subprocess.PIPE)
    result,error = p.communicate()
    sys.stdout.write("Done!\n")
    sys.stdout.write(result.decode('ascii'))

def main():
    (dna_bam, rna_bam, psn, msn, site, analysis_id, token) = get_args()

    sys.stdout.write("Validating DNA and RNA BAM files...\n")
    (new_dna_bam, dna_run_id) = validate_bams(msn, dna_bam, 'dna')
    (new_rna_bam, rna_run_id) = validate_bams(msn, rna_bam, 'rna')

    # Generate a new working directory for our data and move the BAMs there.
    work_dir = gen_wd(psn, msn)
    shutil.move(new_dna_bam, work_dir)
    shutil.move(new_dna_bam + '.bai', work_dir)
    shutil.move(new_rna_bam, work_dir)
    shutil.move(new_rna_bam + '.bai', work_dir)
    os.chdir(work_dir)

    # get IR data
    if not analysis_id:
        if dna_run_id == rna_run_id:
            analysis_id = dna_run_id
        else:
            sys.stderr.write("ERROR: DNA and RNA Analysis IDs do not agree!\n")
            sys.stderr.write("\tDNA: %s\n\tRNA: %s" % (dna_run_id, rna_run_id))
            sys.exit(1)
    get_ir_data(analysis_id, token)

    # Get the VCF file for processing.
    files = os.listdir('vcfs')
    if files:
        if fnmatch.fnmatch(files[0], '*.vcf'):
            vcf = files[0]

    # Generate a MOI report and store it.
    gen_moi_report(msn,'vcfs/' + vcf)

    sys.stdout.write('\nReport generation complete.  Results can be found in {}.\n'.format(os.getcwd()))

if __name__ == '__main__':
    main()
