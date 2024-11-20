#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o errexit

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

options_array=(
  input_directory
  output_directory
  min_sequencing_depth
  max_sequencing_depth
  max_vaf
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'gather_mutect_filter' -- "$@")
eval set -- "${arguments}" while true; do
    case "${1}" in
        --input_directory )
            input_directory="${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --min_sequencing_depth )
            min_sequencing_depth="${2}"; shift 2 ;;
        --max_sequencing_depth )
            max_sequencing_depth="${2}"; shift 2 ;;
        --max_vaf )
            max_vaf="${2}"; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

input_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_5k_filtered"
output_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_5k_filtered_full"
min_sequencing_depth="25"
max_sequencing_depth="100"
max_vaf="0.35"

vcf_list=${input_directory}/vcf_list
find -L "${input_directory}" -type f  `#list all files in ${vcf_directory}` | grep ".*.vcf.gz$" | \
    sort -u  `#sort and remove duplicate names` > ${vcf_list}
vcf_array_length=$(wc -l < ${vcf_list}) #get the number of FASTQs

mkdir -p "${output_directory}/logs"
code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))
sbatch --output "${output_directory}/logs/%A_%a.log" \
    --error "${output_directory}/logs/%A_%a.log" \
    --array "1-${vcf_array_length}" \
    --time 1:00:00 \
    --account sjaiswal \
    --partition batch \
    --cpus-per-task 1 \
    --mem 32G \
    --constraint="nvme" \
    --job-name bcftools_filter \
    "${code_directory}/mutect_filter.sh" \
        --array_file "${vcf_list}" \
        --output_directory ${output_directory} \
        --min_sequencing_depth ${min_sequencing_depth} \
        --max_sequencing_depth ${max_sequencing_depth} \
        --max_vaf ${max_vaf}
