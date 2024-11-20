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
  targets_file
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

input_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_aggregated_chrs"
output_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_all_singletons"
mdust_file="/oak/stanford/groups/sjaiswal/references/mdust_superdups_merged_sorted.bed"
giab_file="/oak/stanford/groups/sjaiswal/references/GRCh38_notinalldifficultregions.bed"
code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))

vcf_list=${input_directory}/vcf_dirs_list
# find -L "${input_directory}" -type f -name "*.vcf.gz" `#list all files in ${vcf_directory}` | \
#     sort -u  `#sort and remove duplicate names` > ${vcf_list}
find -L "${input_directory}" -maxdepth 1 -mindepth 1 -type d  `#list all files in ${vcf_directory}` | grep -v "logs" | \
    sort -u  `#sort and remove duplicate names` > ${vcf_list}


submit_job() {
    # allele_frequency=$1
    # binom_p=$1
    filter_string=$1
    output_name=$2
    output_directory=$3
    vcf_list=$4
    mdust_file=$5
    giab_file=$6
    code_directory=$7

    # output_directory_depth="${output_directory}_af${allele_frequency}"
    # output_directory_depth="${output_directory}_binom_p_${binom_p}"
    # filter_string="AF<=${allele_frequency}"
    output_directory_final="${output_directory}_${output_name}"
    vcf_array_length=$(wc -l < ${vcf_list}) #get the number of FASTQs

    mkdir -p "${output_directory_final}/logs"
    sbatch --output "${output_directory_final}/logs/%A_%a.log" \
        --error "${output_directory_final}/logs/%A_%a.log" \
        --array "1-${vcf_array_length}" \
        --time 6:00:00 \
        --account sjaiswal \
        --partition batch \
        --cpus-per-task 1 \
        --mem 32G \
        --constraint="nvme" \
        --job-name bcftools_concat \
        "${code_directory}/bcftools_concat.sh" \
            --array_file "${vcf_list}" \
            --output_directory ${output_directory_final} \
            --mdust_file "${mdust_file}" \
            --giab_file "${giab_file}" \
            --filter_string "${filter_string}"
}
export -f submit_job

# ad_array=(0.25 0.30 0.35)
# binom_p_array=(0.05 0.1)
# echo ${binom_p_array[@]} | tr ' ' '\n' | xargs -I % bash -c "submit_job % ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}"

submit_job "AD[*:1]>=4 & AF<=0.35" "ad4_af0.35" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}
submit_job "AD[*:1]>=4 & AF<=0.4" "ad4_af0.4" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}
submit_job "AD[*:1]>=4 & AF<=0.45" "ad4_af0.45" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}
submit_job "AD[*:1]>=4 & AF<=0.5" "ad4_af0.5" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}
submit_job "AD[*:1]>=4 & binom(AD[*:1],AD[*:0])<0.05 & AF<=0.5" "ad4_binom_p_0.05" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}
submit_job "AD[*:1]>=4 & binom(AD[*:1],AD[*:0])<0.1 & AF<=0.5" "ad4_binom_p_0.1" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}

# submit_job "AF<=0.35" "af0.35" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}
# submit_job "DP>=25 & AF<=0.35" "dp25_af0.35" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}
# submit_job "AD>=4 & DP>=25 & AF<=0.35" "ad4_dp25_af0.35" ${output_directory} ${vcf_list} ${mdust_file} ${giab_file} ${code_directory}

submit_job "AD[*:1]>=4 & AF<=0.35" "ad4_af0.35" ${output_directory} ${vcf_list} true ${mdust_file} none ${code_directory}
submit_job "AD[*:1]>=4 & AF<=0.35" "ad4_af0.35" ${output_directory} ${vcf_list} false none ${giab_file} ${code_directory}
submit_job "AD[*:1]>=4 & AF<=0.35" "ad4_af0.35" ${output_directory} ${vcf_list} true none ${giab_file} ${code_directory}
submit_job "AD[*:1]>=4 & AF<=0.35" "ad4_af0.35" ${output_directory} ${vcf_list} false ${mdust_file} ${giab_file} ${code_directory}
submit_job "AD[*:1]>=4 & AF<=0.35" "ad4_af0.35" ${output_directory} ${vcf_list} true ${mdust_file} ${giab_file} ${code_directory}
