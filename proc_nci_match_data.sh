#!/usr/bin/env bash
# Pipeline to process NCI-MATCH specimens run by NCI only.  Requires a direct backend connection to MATCHBox and the
# get_mb_data.py script plus the whole variant_review.py package and associated scripts.
# 
# 12/1/2016 - D Sims
#######################################################################################################################
set -e
VERSION='3.4.0_051217'
cwd=$(pwd)

function usage() {
cat << EOT
USAGE: $(basename $0) [-s <sever>] [-a] <psn_list>
    -d    Location of data to download. Options are 'local', 'cloud', and 'matchbox'. DEFAULT: cloud.
    -b    Bucket name needed when downloading data from the cloud instance. DEFAULT: adultmatch.
    -a    Analyze only.  Useful when we already have a directory of PSN_MSN data from get_mb_data.py. This option
          requires no PSN list like when you want to download also (i.e. conventional usage).
    -v    Print version info and exit.
    -h    Print this help text and exit.
EOT
}

while getopts ":had:b:" opt; do
    case $opt in 
        d)
            data_source=$OPTARG
            ;;
        a)
            analyze_only=1
            ;; 
        b)
            bucket=$OPTARG
            ;;
        h)
            #echo "USAGE: $0 [-s <server] [-a] <psn_list>"
            usage
            exit
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

patient_list=$1

# Set default data server
if [[ -z $data_source ]]; then
    #data_source='local'
    data_source='cloud'
fi
if [[ -z $bucket ]]; then 
    bucket='adultmatch'
fi

function check_dir() {
    # Needed just in case we have old analyses in the same dir.  Causes big problems!
    if [[ $(find . -maxdepth 1 -name "PSN*") ]]; then 
        echo "Warning: There are PSN directories in the current working directory already.  Results in those will be 
        purged and regenerated, if you continue."
        read -r -p "Continue [y|N]: " response
        response=${response,,}
        if [[ $response =~ ^(no|n)$ ]]; then
            echo "Exiting.  Move data to safety and try again."
            exit
        elif [[ $response =~ ^(yes|y)$ ]]; then
            echo "Continuing on with fresh analysis of all data."
            rm -rf PSN* MSN*
        else
            echo "Invalid response.  Exiting to be safe"
            exit 1
        fi
    fi
}

function get_data() {
    # Read list, get data
    patient_list=$1
    if [[ -e $patient_list ]]; then
        if [[ $data_source -eq 'cloud' ]]; then
            #get_mb_data.py -s $data_source -a adultmatch -B -b $patient_list 
            get_mb_data.py -d $data_source -a $bucket -s nci -B -b $patient_list 
        else
            get_mb_data.py -s $data_source -B -b $patient_list 
        fi
    else 
        echo "ERROR: no patient list found.  Need to input a file of patient IDs unless running in analyze only mode"
    fi
}

echo "###################################################################"
echo "####    NCI-MATCH Batch Analyzer for NCI Data v$VERSION    ####"
echo "###################################################################"
echo

echo INFO: data source is $data_source
if [[ $analyze_only ]]; then
    echo INFO: Running analysis only on pre-downloaded data.
else
    check_dir && get_data $patient_list
fi

# Get list of PSN_MSNs 
dir_regex='^./PSN[0-9]+_MSN[0-9].*'
echo -n "Generating a list of MSNs from PSNs..."
samples=( $(find . -maxdepth 1 -name "PSN*" | sed 's/\.\///') )
echo "Done!"

# Collect the BAM files, scrap the rest.
echo -n "Moving BAM files to current directory to stage them for analysis..."
find -E . -maxdepth 2 -regex $dir_regex -name "*bam" -exec mv {} $cwd \; && rm -r PSN*
echo "Done!"

# from the data we just generated and run variant_review on them.
count=${#samples[@]}
echo "Running variant_review on data ($count samples to process)..."
iter=1
trap 'exit' INT
for s in ${samples[@]}; do
    echo -e "\t[$iter/$count]  Processing $s."
    patient=(${s//_/ })
    variant_review.py -m ${patient[1]} -p ${patient[0]} ${patient[1]}*[dr]* 2>&1 > /dev/null
    iter=$((iter+1))
done
echo "Done!"
echo "Results are ready for review."
