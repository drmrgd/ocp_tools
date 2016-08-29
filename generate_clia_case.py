#!/usr/bin/python
# Auto generate an entry in the clia_case_list.csv file in order to make adding the batchwise MSNs a little easier.
import sys
import os
import re
import subprocess

version = '0.0.5_081116'

def main():
    if len(sys.argv) < 3:
        sys.stderr.write('ERROR: not enough args passed to script!\n')
        sys.stderr.write('USAGE: {} <run_directory> <clia_case_list.csv>\n'.format(sys.argv[0]))
        sys.exit(1)
    rundir,caselist = sys.argv[1:]
    run_name = re.match('^Auto_user_(.*?\w{3})_\d+_\d+\/?$',rundir).group(1)
    msn_list = get_msn_list(rundir)
    next_casenum = gen_casenum(caselist)
    gen_case_list_string(run_name,msn_list,next_casenum,caselist)
    return

def gen_casenum(caselist):
    with open(caselist) as fh:
        last_entry = fh.readlines()[-1]
        data = last_entry.split(',')
        (project,number) = data[0].split('-')
        next_num = "{0:05d}".format(int(number) + 1)
        next_case = project +'-'+ str(next_num)
        return next_case

def get_msn_list(rundir):
    params_file = rundir + '/ion_params_00.json'
    p = subprocess.Popen(['sampleKeyGen.pl', '-p', '-r', '-f', params_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result,error = p.communicate()
    msn_list = dict([x.split('\t') for x in result.split('\n') if x]) 
    uniq_msn_set = set(v for v in msn_list.values())
    return [x for x in uniq_msn_set]

def gen_case_list_string(run_name,msns,new_casenum,caselist):
    msn_list = ','.join(msns)
    new_entry=','.join([new_casenum,msn_list,run_name])
    print "Adding new entry to case list file:\n{}".format(new_entry)
    with open(caselist,'a') as fh:
        fh.write(','.join([new_casenum,msn_list,run_name]))
        fh.write('\n')

if __name__ == '__main__':
    main()
