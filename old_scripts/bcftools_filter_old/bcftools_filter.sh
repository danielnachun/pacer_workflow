#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

check_for_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ${file_path} != "none" ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
options_array=(
    array_file
    output_directory
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bcftools_filter' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --array_file )
            array_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

line_number=${SLURM_ARRAY_TASK_ID}
vcf_list="$(sed "${line_number}q; d" "${array_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

chr_name=$(bcftools index -s $(head -n 1 ${vcf_list}) | cut -f 1)
dir_name=$(basename $(dirname ${vcf_list}))
mkdir -p ${output_directory}/${dir_name}

# Merge all samples for a given chromosome
# Add AC tag to get allele count
# Only keep biallelic variants which pass all filtering conditions
bcftools merge -Ou -l ${vcf_list} | \
    bcftools +fill-tags -Ou -- --tags AC | \
    bcftools view --min-alleles 2 --max-alleles 2 \
        --include "AD[*:1]>=4 & AF<=0.4 & DP<=100 & FILTER='PASS' & INFO/AC=1" \
        -o ${output_directory}/${dir_name}/${chr_name}.filter.vcf.gz
