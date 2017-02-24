#!/usr/bin/python
import sys
import os
import json
import argparse
from collections import defaultdict
from pprint import pprint as pp

version = '0.9.0_022417'

class FusionPanel(object):
    def __init__(self,fusions_json):
        self.fusions_json = fusions_json

    def __getitem__(self,key):
        return self.fusion_panel[key]

    def __iter__(self):
        return self.fusion_panel.itervalues()

    @classmethod
    def read_json(cls,fusion_json):
        '''Read in a config file of params to use in this program'''
        with open(fusion_json) as fh:
            data = json.load(fh)
        return data


def get_args():
    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 100, width=200),
        description='''
        <Program Description>
        ''',
        version = '%(prog)s  - ' + version,
        )
    parser.add_argument('vcf', metavar='<VCF(s)>', nargs='+', help='OCAv3 VCF file to parse.')
    parser.add_argument('-o', '--output', default='output.txt', metavar='<outfile>', help='Custom output file (DEFAULT: %(default)s)')
    args = parser.parse_args()
    return args

def parse_line(line,assay_type,panel):
    elems = line.split(';')
    for e in elems:
        if e.startswith('ID'):
            gene=get_val(e)
        elif e.startswith('READ_COUNT'):
            counts=get_val(e)
    return gene,counts,panel[assay_type][gene]

def get_val(string):
    return string.split('=')[1]

def trim_string(line):
    elems=line.split()
    return 'ID={};{}'.format(elems[2],elems[7])

def parse_control_results(data,pool):
    return sum([int(data[x][0]) for x in data if data[x][1] == pool])

def read_vcf(vcf,panel):
    '''Read in VCF and return expr control and gene expression results as a dict of lists for downstream parsing'''
    mapped_reads = ''
    fusion_calls = defaultdict(dict)

    with open(vcf) as fh:
        for line in fh:
            if line.startswith('##TotalMappedFusionPanelReads'):
                mapped_reads = line.split('=')[1].rstrip()
            if line.startswith('#CHROM'):
                sample_name = line.split()[-1]
            elif 'SVTYPE=ExprControl' in line:
                gene,counts,pool = parse_line(trim_string(line),'ExpressionControl',panel)
                fusion_calls['expr_ctrl'][gene] = (counts,pool)
            elif 'SVTYPE=GeneExpression' in line:
                gene,counts,pool = parse_line(trim_string(line),'GeneExpression',panel)
                fusion_calls['gene_expr'][gene] = (counts,pool)
    return fusion_calls,mapped_reads,sample_name

def print_data(data):
    print(','.join(['sample_name','mapped_reads','pool1_ec_reads','pool2_ec_reads','pool1_ge_reads','pool2_ge_reads','pool1_total','pool2_total']))
    for sample in sorted(data):
        sum1 = sum(map(int,[data[sample][1],data[sample][3]]))
        sum2 = sum(map(int,[data[sample][2],data[sample][4]]))
        new_list = [x for x in [sample]+data[sample]+[sum1,sum2]]
        print(','.join(str(x) for x in new_list))

def raw_print(data):
    print(','.join(['gene','counts','pool']))
    for assay in data:
        for gene in data[assay]:
            print(','.join(map(str,[gene,data[assay][gene][0],data[assay][gene][1]])))
                
def main():
    args = get_args()
    fusion_panel_json = os.path.dirname(os.path.realpath(__file__)) + '/fusion_panel.json'
    fusion_panel = FusionPanel.read_json(fusion_panel_json)
    
    results = {}
    for vcf in args.vcf:
        control_data,mapped_reads,sample_name = read_vcf(vcf,fusion_panel)

        pool1_ec_sum = parse_control_results(control_data['expr_ctrl'],'pool1')
        pool2_ec_sum = parse_control_results(control_data['expr_ctrl'],'pool2')
        pool1_ge_sum = parse_control_results(control_data['gene_expr'],'pool1')
        pool2_ge_sum = parse_control_results(control_data['gene_expr'],'pool2')
        results[sample_name] = [mapped_reads, pool1_ec_sum, pool2_ec_sum, pool1_ge_sum, pool2_ge_sum]

        raw_print(control_data)

    # print("Total mapped reads: {}".format(mapped_reads))
    # print("pool1 expression control sum: {}\npool2 exppression conrol sum: {}\n".format(pool1_ec_sum,pool2_ec_sum))
    # print("pool1 gene expression sum: {}\npool2 gene expression sum: {}\n".format(pool1_ge_sum,pool2_ge_sum))

    print_data(results)
    
if __name__ == '__main__':
    main()
