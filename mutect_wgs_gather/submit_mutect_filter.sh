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
  bravo_variants
  reference_genome
  sequence_dictionary
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'gather_mutect_filter' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --input_directory )
            input_directory="${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --bravo_variants )
            bravo_variants="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sequence_dictionary )
            sequence_dictionary="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

reference_genome="/oak/stanford/groups/sjaiswal/references/cloud_references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
sequence_dictionary="/oak/stanford/groups/sjaiswal/references/cloud_references/GRCh38_full_analysis_set_plus_decoy_hla.dict"
bravo_variants="/oak/stanford/groups/sjaiswal/references/BRAVO_TOPMed_Freeze_8.vcf.gz"
input_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_next_batch2/outputs_full_intervals"
output_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_next_batch_filtered"

vcf_list=${input_directory}/vcf_list
find -L "${input_directory}" -maxdepth 1 -mindepth 1 -type d  `#list all files in ${vcf_directory}` | \
    sort -u  `#sort and remove duplicate names` > ${vcf_list}
vcf_array_length=$(wc -l < ${vcf_list}) #get the number of FASTQs

mkdir -p "${output_directory}/logs"
code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))
sbatch --output "${output_directory}/logs/%A_%a.log" \
    --error "${output_directory}/logs/%A_%a.log" \
    --array "1-${vcf_array_length}" \
    --time 12:00:00 \
    --account sjaiswal \
    --partition batch \
    --cpus-per-task 1 \
    --mem 32G \
    --constraint="nvme" \
    --job-name gather_mutect_filter \
    "${code_directory}/mutect_filter.sh" \
        --array_file "${vcf_list}" \
        --output_directory ${output_directory} \
        --bravo_variants ${bravo_variants} \
        --reference_genome ${reference_genome} \
        --sequence_dictionary ${sequence_dictionary}
