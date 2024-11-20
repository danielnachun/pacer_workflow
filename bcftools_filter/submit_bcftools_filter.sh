#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o errexit

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
  vcf_list_directory
  output_directory
  chunk_size
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'submit_bcftools_merge' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --input_directory )
            input_directory="${2}"; shift 2 ;;
        --vcf_list_directory )
            vcf_list_directory="${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

input_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_aggregated_all"
vcf_list_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_aggregated_chrs_vcf_lists"
output_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_aggregated_chrs"

mkdir -p ${vcf_list_directory}

chrs=$(bcftools index -s $(find -L ${input_directory} -type f -name '*.vcf.gz' | grep -v 'chr' | sort -u | head -n 1) | cut -f 1)
echo ${chrs} | tr ' ' '\n' | xargs --replace=% bash -c "find -L ${input_directory} -type f -name '*%.vcf.gz' | sort -u > ${vcf_list_directory}/%_vcfs"

vcf_lists=$(dirname ${input_directory})/vcf_lists
find -L "${vcf_list_directory}" -type f `#list all files in ${vcf_list_directory}` | \
    sort -u > ${vcf_lists}
vcf_array_length=$(wc -l < ${vcf_lists}) #get the number of VCF lists

mkdir -p "${output_directory}/logs"
code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))
sbatch --output "${output_directory}/logs/%A_%a.log" \
    --error "${output_directory}/logs/%A_%a.log" \
    --array "1-${vcf_array_length}" \
    --time 6:00:00 \
    --account sjaiswal \
    --partition batch \
    --cpus-per-task 1 \
    --mem 32G \
    --constraint="nvme" \
    --job-name bcftools_filter \
    "${code_directory}/bcftools_filter.sh" \
        --array_file "${vcf_lists}" \
        --output_directory ${output_directory} 
