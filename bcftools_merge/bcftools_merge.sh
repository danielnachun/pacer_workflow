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
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bcftools_merge' -- "$@")
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

# Merge 100 samples at a time into single VCFs
# Then extract each chromosome and only keep INFO and FORMAT tags needed later
bcftools merge --force-samples -l ${vcf_list} -o ${output_directory}/aggregated_vcf_${line_number}.vcf.gz
tabix -p vcf ${output_directory}/aggregated_vcf_${line_number}.vcf.gz
bcftools index -s ${output_directory}/aggregated_vcf_${line_number}.vcf.gz | \
    cut -f 1 | \
    xargs -I % bash -c "bcftools view ${output_directory}/aggregated_vcf_${line_number}.vcf.gz % | \
        bcftools annotate -x 'INFO,^FORMAT/GT,FORMAT/DP,FORMAT/AF,FORMAT/AD' \
        -o ${output_directory}/aggregated_vcf_${line_number}.%.vcf.gz; \
        tabix -p vcf ${output_directory}/aggregated_vcf_${line_number}.%.vcf.gz"
