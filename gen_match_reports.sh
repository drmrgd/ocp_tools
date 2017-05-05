#!/bin/bash
# Generate all the report files that I want to send out for the MATCH Validation
# study.
# 
# 3/13/15 - D Sims
###################################################################################
version='2.0.0_050517'

if [[ $1 == '-h' ]]; then
    echo "USAGE: gen_match_reports.sh"
    echo "Run this in a directory of VCFs and produce a series of NCI-MATCH reports for analysis.  No args to pass."
    exit
fi

# Generate MOI reports
echo "Generating MATCH MOI Reports..."
if [[ ! -d 'moi_reports' ]]; then
    mkdir ./moi_reports
fi

parallel "match_moi_report.pl -c7 -r1000 -o {= s/-DNA.*/_moi_report.txt/ =} 2>&1 > /dev/null {}" ::: *vcf
mv *moi_report.txt ./moi_reports/
collate_moi_reports.py -o 'collated_moi_reports.csv' --cn 7 --reads 1000 *vcf
echo "Done with MOI Reports!"

# Generate VCF metrics data
echo "Generating QC Report from VCF files..."
get_metrics_from_vcf.py -o vcf_metrics.txt *vcf
echo 'Done!'
