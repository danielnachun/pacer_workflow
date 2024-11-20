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
  mdust_file
  giab_file
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'gather_mutect_filter' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --input_directory )
            input_directory="${2}"; check_for_directory ${1} ${2}; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --mdust_file )
            mdust_file="${2}"; check_for_file ${1} ${2}; shift 2 ;;
        --giab_file )
            giabl_file="${2}"; check_for_file ${1} ${2}; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

input_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_aggregated_chrs"
output_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_all_singletons"
mdust_file="/oak/stanford/groups/sjaiswal/references/mdust_superdups_merged_sorted.bed"
giab_file="/oak/stanford/groups/sjaiswal/references/GRCh38_notinalldifficultregions.bed"

code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))

mkdir -p "${output_directory}/logs"
sbatch --output "${output_directory}/logs/%A_%a.log" \
    --error "${output_directory}/logs/%A_%a.log" \
    --array "1-${vcf_array_length}" \
    --time 6:00:00 \
    --account sjaiswal \
    --partition batch \
    --cpus-per-task 1 \
    --mem 32G \
    --constraint="nvme" \
    --job-name bcftools_concat \
    "${code_directory}/bcftools_concat.sh" \
        --input_directory "${input_directory}" \
        --output_directory ${output_directory} \
        --mdust_file "${mdust_file}" \
        --giab_file "${giab_file}"
