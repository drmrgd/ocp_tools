#!/usr/bin/env python3
# Read in a list of MSNs, and get + delink the data for use in other experiments.  Maintain original
# link list until everything looks OK, and then manually remove the original data.
#
# 11/10/2016 - D Sims
########################################################################################################
import sys
import re
import os
import random
import argparse
import shutil
import datetime
import subprocess
from collections import defaultdict
from pprint import pprint as pp

version = '1.3.0_030817'
cwd = os.getcwd()

def get_args():
    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 100, width=200),
        description=
        '''
        Program to delink MATCH Data.  Based on input in the form of an MSN list, plus a directory containing the 
        corresponding BAM and VCF files (usually derived from 'get_mb_data.py'), delink an run data, and return a 
        sample key file containing the original and new identifiers (for the temporary purpose of renaming tubes or 
        whatever), and a new dataset with relevent identifiers removed from the file.  
        ''',
        )
    parser.add_argument('sample_file', metavar='<msn_list.file>', 
            help='Flat file list of MSNs corresponding with the samples you wish to delink')
    parser.add_argument('sample_dirs', metavar='<sample_directories>', nargs='+', 
            help='Directories that contain the BAM and VCF files matching the sample list (e.g. PSN12345_MSN6789)')
    parser.add_argument('-p','--prefix', metavar='<string>', default='Sample',
            help='Prefix for new sample name that would preceed the randomized number string (DEFAULT: "%(default)s")')
    parser.add_argument('-v', '--version', action='version', version = '%(prog)s - ' + version) 
    args = parser.parse_args()
    return args

def read_sample_list(sample_file,prefix):
    '''Create a randomized sample name and a sampleKey.txt file for temporary renaming purposes'''
    with open(sample_file) as fh:
        samples = ['MSN'+i.rstrip('\n').lstrip('MSN') for i in fh] # make sure we have a MSN designator in string.
        
    rnd = gen_rand_name(len(samples))
    name = [prefix +'-'+ str(x).zfill(4) for x in rnd]
    new_names = dict(zip(samples,name))

    with open('sampleKey.txt', 'w') as outfh:
        for i in new_names:
            outfh.write('{},{}\n'.format(i,new_names[i]))
    return new_names 

def gen_rand_name(num_elems):
    return random.sample(range(0,9999),num_elems)

def check_manifest(dirs,samples):
    '''Based on the directories loaded and the sample list, make sure all match'''
    orig_samples = samples.keys()
    count = defaultdict(int)
    skipped_dirs = []
    
    print('Validating manifest of samples and directories to process...',end='')
    for d in dirs:
        msn = d.rstrip('/').split('_')[1]
        if msn in orig_samples:
            count[msn] += 1
        else:
            sys.stderr.write("WARN: {} does not have an entry in the original sample manifest file.  Skipping this dataset (NOTE: All skipped runs found in 'skipped_dirs' directory).\n".format(d))
            skipped_dirs.append(d)
            del dirs[d]

    if skipped_dirs:
        os.mkdir('skipped_expts')
        for d in skipped_dirs:
            new_location = os.path.join('skipped_expts',d)
            shutil.move(d,new_location)

    extra_samples = []
    for s in samples:
        if not count[s]:
            sys.stderr.write("WARN: {} is in the sample list, but the data was not found! Removing this one from list and skipping.\n".format(s))
            extra_samples.append(s)

    if extra_samples:
        samples = {s : samples[s] for s in samples if s not in extra_samples}
    print('Done!')
    return samples, dirs

def get_header(line):
    return line.split('=')[0]

def time():
    now = str(datetime.datetime.now().strftime('%Y%m%d'))
    utc = str(datetime.datetime.utcnow())
    return now,utc.replace(' ', 'T')

def cleanup(sample_id):
    wd = os.getcwd()
    all_files = os.listdir(wd)
    for f in all_files:
        if not f.startswith(sample_id):
            os.remove(f)

def proc_vcf(vcf,delinked_id):
    '''Read in VCF file and substitute the appropriate lines with new data to delink the specimen'''
    now,utc = time()
    sys.stdout.write('\tDelinking VCF file %s...' % vcf)
    new_vcf = delinked_id + '.vcf'
    out_fh = open(new_vcf, 'w')
    
    with open(vcf) as vcf_fh:
        for line in vcf_fh:
            if line.startswith('##fileDate'):
                new_line = get_header(line) + '=' +  now + ' (delinked)'
                out_fh.write("{}\n".format(new_line))
            elif line.startswith('##fileUTCtime'):
                new_line = get_header(line) +'='+ utc + ' (delinked)'
                out_fh.write("{}\n".format(new_line))
            elif line.startswith('#CHROM'):
                elems = line.split()
                elems[9] = delinked_id
                out_fh.write('\t'.join(elems))
                out_fh.write('\n')
            else:
                out_fh.write(line)
    out_fh.close()
    sys.stdout.write('Done!\n')

def proc_bam(bam,orig_id,delinked_id):
    '''Read BAM header in and edit appropriate information to delink specimen'''
    if bam.endswith('rna.bam'):
        new_bam = delinked_id + '_rna.bam'
    elif bam.endswith('dna.bam'):
        new_bam = delinked_id + '_dna.bam'

    now,utc = time()
    bam_fh = open(new_bam, 'w')

    sys.stdout.write('\tDelinking and Reheadering BAM file: {}...'.format(bam))
    sys.stdout.flush()
    try:
        p = subprocess.check_output(['samtools', 'view', '-H', bam], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as err:
        print('ERROR: Can not read BAM header: {} ({})!'.format(err,err.returncode))
    header = p.decode('ascii').split('\n')

    # Create a temporary header that we can use to fix BAM.
    with open('tmp.header', 'a') as out_fh:
        for line in header:
            line = line.replace(orig_id,delinked_id)
            line = re.sub('DT:.*?\t','DT:%s\t' % utc,line)
            line = re.sub('CL:.*?$','CL:<delinked>',line)
            out_fh.write("%s\n" % line)

    # Re-header the BAM file
    try:
        subprocess.run(['samtools', 'reheader', '-P', 'tmp.header', bam], stdout=bam_fh, check=True)
        os.remove('tmp.header')
    except subprocess.CalledProcessError as err:
        sys.stderr.write("ERROR: failed to re-header the BAM file: {}!\n".format(err,err.returncode))
        sys.exit(1)
    sys.stdout.write('Done!\n')

def delink_data(sample_list,dirs):
    '''For each elem in the sample list dict, find appropriate dict, read in VCF file and change, read in BAM file 
       change. Move all original data to a copies dir to make sure we have what we need before we finish'''

    # Create a place to store original data before we expunge it.
    if not os.path.exists('orig_data'):
        os.mkdir('orig_data')
    
    count = 0
    for sample in sample_list:
        count += 1
        print('  [{}/{}] Processing sample: {}'.format(count,len(sample_list),sample))
        for d in dirs:
            if d.endswith(sample):
                if not os.path.exists('orig_data/' + d):
                    sys.stdout.write('\tCopying {} to backup location...'.format(d))
                    sys.stdout.flush()
                    shutil.copytree(d,'orig_data/'+d)
                    sys.stdout.write('Done!\n')

                files = os.listdir(d)
                os.chdir(d)
                for f in files:
                    if f.startswith(sample) and f.endswith('vcf'):
                        proc_vcf(f,sample_list[sample])
                    elif f.startswith(sample) and f.endswith('bam'):
                        proc_bam(f,sample,sample_list[sample])
                cleanup(sample_list[sample])
                os.chdir(cwd)
                os.rename(d,sample_list[sample])
                print()
                continue

if __name__=='__main__':
    args = get_args()

    sample_list = read_sample_list(args.sample_file,args.prefix)

    dirs = [d.replace('/','') for d in args.sample_dirs]
    final_samplelist,final_dirlist = check_manifest(dirs,sample_list)

    sys.stdout.write("Delinking {} files based on input list.\n".format(len(final_samplelist)))
    sys.stdout.flush()
    delink_data(final_samplelist, final_dirlist)
    sys.stdout.write("All MATCH data for manifest is now delinked. Original data stored in 'orig_data' dir and must be deleted manually\n")
