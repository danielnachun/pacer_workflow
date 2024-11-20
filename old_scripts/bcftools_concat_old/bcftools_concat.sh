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
    mdust_file
    giab_file
    filter_string
    cttc
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'somatic_variant_pipeline' -- "$@")
eval set -- "${arguments}"

giab_file="none"
mdust_file="none"
cttc="false"

while true; do
    case "${1}" in
        --array_file )
            array_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --mdust_file )
            mdust_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --giab_file )
            giab_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --filter_string )
            filter_string="${2}"; shift 2 ;;
        --cttc )
            cttc="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

line_number=${SLURM_ARRAY_TASK_ID}
vcf_dir="$(sed "${line_number}q; d" "${array_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

#vcf_file="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_5k_aggregated_saturation/18/chr1.filter.vcf.gz"
#output_directory="/oak/stanford/groups/sjaiswal/dnachun/ukbb/chip_carriers_5k_singletons_test"
#targets_file="/oak/stanford/groups/sjaiswal/references/mdust_superdups_merged_sorted.bed"

# dir_name=$(basename $(dirname ${vcf_dir}))
dir_name=$(basename ${vcf_dir})
mkdir -p ${output_directory}/${dir_name}
# chr_name=$(echo $(basename ${vcf_file}) | sed -e 's/\..*$//g')
# tabix -f -p vcf ${vcf_file}
# bcftools view ${vcf_file} --include "${filter_string}" --output ${output_directory}/${dir_name}/${chr_name}_singletons.vcf.gz
# tabix -p vcf ${output_directory}/${dir_name}/${chr_name}_singletons.vcf.gz
bcftools concat -Ou --file-list <(find ${vcf_dir} -name '*.vcf.gz' | grep -v -E "chrX|chrY") | bcftools view --include "${filter_string}" -o ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz
tabix -p vcf ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz

bcftools stats --samples - ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz > ${output_directory}/${dir_name}/aggregated_stats.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.txt

# bcfview_args=()
# bcftools_cmds=()
# bcfview_args=""
# bcftools_cmds=""
# stats_file_name="aggregated_stats"
# singletons_file_name="aggregated_singletons"

# if [[ ${cttc} == "true" ]]; then
    # bcfview_include="(REF='C' & ALT='T') | (REF='T' & ALT='C') | (REF='A' & ALT='G') | (REF='G' & ALT='A')"
    # bcfview_args+=(--include "${bcfview_include[@]}")
#     stats_file_name="${stats_file_name}.cttc"
#     singletons_file_name="${singletons_file_name}.cttc"
# fi

# if [[ ${mdust_file} != "none" ]]; then
#     bcftools_cmds="${bcftools_cmds} | bcftools view -Ou --targets-file '^${mdust_file}'"
#     stats_file_name="${stats_file_name}.mdust"
#     singletons_file_name="${singletons_file_name}.mdust"
# fi

# if [[ ${giab_file} != "none" ]]; then
#     bcftools_cmds="${bcftools_cmds} | bcftools view -Ou --targets-file '${giab_file}'"
#     stats_file_name="${stats_file_name}.giab"
#     singletons_file_name="${singletons_file_name}.giab"
# fi

# bcfview_cmd="bcftools view -Ou ${bcfview_args} ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz"
# ${bcfview_cmd} ${bcftools_cmds} | \
#     bcftools stats --samples - > ${output_directory}/${dir_name}/${stats_file_name}.txt
# grep "PSC" < ${output_directory}/${dir_name}/${stats_file_name}.txt | \
#     tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/${singletons_file_name}.txt

bcftools view -Ou --include '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="A" & ALT="G") | (REF="G" & ALT="A")' \
    ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
    bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.cttc.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.cttc.txt | \
    tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.cttc.txt

bcftools view -Ou --targets-file "^${mdust_file}" ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
    bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.mdust.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.mdust.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.mdust.txt

bcftools view -Ou --include '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="A" & ALT="G") | (REF="G" & ALT="A")' ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
    bcftools view -Ou --targets-file "^${mdust_file}" | \
    bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.mdust.cttc.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.mdust.cttc.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.mdust.cttc.txt

# GIAB filtering
bcftools view -Ou --targets-file "${giab_file}" ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
    bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.giab.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.giab.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.giab.txt

bcftools view -Ou --include '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="A" & ALT="G") | (REF="G" & ALT="A")' \
    ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
    bcftools view -Ou --targets-file "${giab_file}" | \
    bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.giab.cttc.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.giab.cttc.txt | \
    tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.giab.cttc.txt

bcftools view -Ou --targets-file "${giab_file}" ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
    bcftools view -Ou --targets-file "^${mdust_file}" | \
    bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.giab.mdust.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.giab.mdust.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.giab.mdust.txt

bcftools view -Ou --include '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="A" & ALT="G") | (REF="G" & ALT="A")' ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
    bcftools view -Ou --targets-file "${giab_file}" | \
    bcftools view -Ou --targets-file "^${mdust_file}" | \
    bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.giab.mdust.cttc.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.giab.mdust.cttc.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.giab.mdust.cttc.txt
