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

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
options_array=(
    array_file
    output_directory
    min_sequencing_depth
    max_sequencing_depth
    max_vaf
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'somatic_variant_pipeline' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --array_file )
            array_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --min_sequencing_depth )
            min_sequencing_depth="${2}"; shift 2 ;;
        --max_sequencing_depth )
            max_sequencing_depth="${2}"; shift 2 ;;
        --max_vaf )
            max_vaf="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

line_number=${SLURM_ARRAY_TASK_ID}
filtered_vcf="$(sed "${line_number}q; d" "${array_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

# Filter variants that do not pass filter
filtered_vcf_filter_pass="${output_directory}/$(basename ${filtered_vcf//.vcf.gz/_full.vcf.gz})"
if [[ ! -f "${filtered_vcf_filter_pass}" ]]; then
    filter_array=(
            'FILTER="PASS"'
            "FORMAT/DP>=${min_sequencing_depth}"
            "FORMAT/DP<=${max_sequencing_depth}"
            "AF<=${max_vaf}"
        )
    filter_string=$(printf '%s & ' "${filter_array[@]}" | sed 's/ & $//g')
    bcftools view "${filtered_vcf}" --include "${filter_string}" --output "${filtered_vcf_filter_pass}"
    tabix --force --preset vcf "${filtered_vcf_filter_pass}"
fi
