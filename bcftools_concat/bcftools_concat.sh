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
    input_directory
    output_directory
    mdust_file
    giab_file
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'somatic_variant_pipeline' -- "$@")
eval set -- "${arguments}"
# Set default values for some arguments
giab_file="none"
mdust_file="none"

while true; do
    case "${1}" in
        --input_directory )
            input_directory="${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --mdust_file )
            mdust_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --giab_file )
            giab_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --filter_string )
            filter_string="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# Concatenate all chromosomes into one VCF
mkdir -p ${output_directory}
bcftools concat -Ou --file-list <(find ${input_directory} -name '*.vcf.gz' | \
    grep -v -E "chrX|chrY") | \
    bcftools view --include "${filter_string}" -o ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz
tabix -p vcf ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz

# Get counts for all remaining singletons
bcftools stats --samples - ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz > ${output_directory}/${dir_name}/aggregated_stats.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.txt

# Get counts for C>T/T>C singletons
bcftools view -Ou --include '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="A" & ALT="G") | (REF="G" & ALT="A")' \
    ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
    bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.cttc.txt
grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.cttc.txt | \
    tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.cttc.txt

# We are not using these filters anymore, but the code is kept here for now
# Get counts for singletons excluding mdust/superDups
# bcftools view -Ou --targets-file "^${mdust_file}" ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
#     bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.mdust.txt
# grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.mdust.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.mdust.txt

# Get counts for C>T/T>C singletons excluding mdust/superDups intervals
# bcftools view -Ou --include '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="A" & ALT="G") | (REF="G" & ALT="A")' ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
#     bcftools view -Ou --targets-file "^${mdust_file}" | \
#     bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.mdust.cttc.txt
# grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.mdust.cttc.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.mdust.cttc.txt

# Get counts for singletons only in "good" regions from Genome in a Bottle
# bcftools view -Ou --targets-file "${giab_file}" ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
#     bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.giab.txt
# grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.giab.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.giab.txt

# Get counts for C>T/T>C singletons only in "good" regions from Genome in a Bottle
# bcftools view -Ou --include '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="A" & ALT="G") | (REF="G" & ALT="A")' \
#     ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
#     bcftools view -Ou --targets-file "${giab_file}" | \
#     bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.giab.cttc.txt
# grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.giab.cttc.txt | \
#     tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.giab.cttc.txt

# Get counts for singletons only in "good" regions from Genome in a Bottle excluding mdust/superDups intervals
# bcftools view -Ou --targets-file "${giab_file}" ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
#     bcftools view -Ou --targets-file "^${mdust_file}" | \
#     bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.giab.mdust.txt
# grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.giab.mdust.txt | tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.giab.mdust.txt

# Get counts for C>T/T>C singletons only in "good" regions from Genome in a Bottle excluding mdust/superDups intervals
# bcftools view -Ou --include '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="A" & ALT="G") | (REF="G" & ALT="A")' ${output_directory}/${dir_name}/aggregated_singletons.vcf.gz | \
#     bcftools view -Ou --targets-file "${giab_file}" | \
#     bcftools view -Ou --targets-file "^${mdust_file}" | \
#     bcftools stats --samples - > ${output_directory}/${dir_name}/aggregated_stats.giab.mdust.cttc.txt
# grep "PSC" < ${output_directory}/${dir_name}/aggregated_stats.giab.mdust.cttc.txt | \
#     tail -n +2 | cut -f 3,11 | sed -e 's/\[[0-9]*\]//g' > ${output_directory}/${dir_name}/aggregated_singletons.giab.mdust.cttc.txt
