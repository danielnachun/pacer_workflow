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
    code_directory
    reference_genome
    sequence_dictionary
    bravo_variants
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
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sequence_dictionary )
            sequence_dictionary="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --bravo_variants )
            bravo_variants="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

line_number=${SLURM_ARRAY_TASK_ID}
original_input_directory="$(sed "${line_number}q; d" "${array_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID
sample_name=$(basename "${original_input_directory}")

final_output_directory="${output_directory}/${sample_name}"
mkdir -p ${final_output_directory}

output_directory=${TMPDIR}/${sample_name}
mkdir -p ${output_directory}
input_directory=${TMPDIR}/${sample_name}
mkdir -p ${input_directory}
rsync -Prhltv ${original_input_directory}/ ${input_directory}
reference_directory=${TMPDIR}/references
mkdir -p ${reference_directory}
rsync -Prhltv $(dirname ${reference_genome})/ ${reference_directory}
rsync -Prhltv ${bravo_variants} ${reference_directory}
rsync -Prhltv ${bravo_variants}.tbi ${reference_directory}

# Gather pileup summaries
if [[ ! -f "${output_directory}/${sample_name}_pileups.table" ]]; then
    pileup_tables=$(find "${input_directory}/pileups" -type f | grep -E ".*_pileups.table$" | sort -V | sed -e 's/^/--I /g' | tr '\n' ' ')
    read -r -a pileup_tables_array <<< "${pileup_tables}"
    gatk GatherPileupSummaries \
        --sequence-dictionary "${reference_directory}/$(basename ${sequence_dictionary})" \
        "${pileup_tables_array[@]}" \
        -O "${output_directory}/${sample_name}_pileups.table"
fi

# Gather VCFs
if [[ ! -f "${output_directory}/${sample_name}_mutect2.vcf" ]]; then
    vcf_files=$(find "${input_directory}/vcfs" -type f | grep -E ".*.vcf$" | sort -V | sed -e 's/^/--INPUT /g' | tr '\n' ' ')
    read -r -a vcf_files_array <<< "${vcf_files}"
    gatk MergeVcfs \
        "${vcf_files_array[@]}" \
        --OUTPUT "${output_directory}/${sample_name}_mutect2.vcf"
fi

# Gather VCF stats
if [[ ! -f "${output_directory}/${sample_name}_mutect2.vcf.stats" ]]; then
    vcf_stats=$(find "${input_directory}/vcfs" -type f | grep -E ".*.vcf.stats$" | sort -V | sed -e 's/^/--stats /g' | tr '\n' ' ')
    read -r -a vcf_stats_array <<< "${vcf_stats}"
    gatk MergeMutectStats \
        "${vcf_stats_array[@]}" \
        --output "${output_directory}/${sample_name}_mutect2.vcf.stats"
fi

# Learn read orientatation bias model
if [[ ! -f "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz" ]]; then
    echo "Learning read orientation bias model..."
    f1r2_files=$(find "${input_directory}/f1r2" -type f | sed -e 's/^/-I /g' | tr '\n' ' ')
    read -r -a f1r2_files_array <<< "${f1r2_files}"
    gatk LearnReadOrientationModel \
        "${f1r2_files_array[@]}" \
        --output "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz"
    echo "...read orientation bias model trained."
fi
    
# Calculate cross sample contamination using pileups at known common germline variants.
if [[ ! -f "${output_directory}/${sample_name}_contamination.table" ]]; then
    echo "Getting contamination rate with CalculateContamination..."
    gatk CalculateContamination \
        --input "${output_directory}/${sample_name}_pileups.table" \
        --output "${output_directory}/${sample_name}_contamination.table"
    echo "...contamination rate calculated."
fi

# Add FILTER column to Mutect2 VCF to identify variants which pass or fail filters.
filtered_vcf="${output_directory}/${sample_name}_mutect2_filtered.vcf"
if [[ ! -f "${filtered_vcf}" ]]; then
    echo "Filtering somatic variants with FilterMutectCalls..."
    gatk FilterMutectCalls \
        --variant "${output_directory}/${sample_name}_mutect2.vcf" \
        --output "${filtered_vcf}" \
        --contamination-table "${output_directory}/${sample_name}_contamination.table" \
        --ob-priors  "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz" \
        --reference "${reference_directory}/$(basename ${reference_genome})"
    echo "...somatic variants filtered."
fi

# Only keep biallelic SNVs
filtered_vcf_filter_pass="${filtered_vcf//.vcf/_filter_pass.vcf.gz}"
if [[ ! -f "${filtered_vcf_filter_pass}" ]]; then
    filter_string='TYPE="snp"'
    bcftools view "${filtered_vcf}" --include "${filter_string}" --min-alleles 2 --max-alleles 2 --output "${filtered_vcf_filter_pass}"
    tabix --force --preset vcf "${filtered_vcf_filter_pass}"
fi

# Remove variants which are BRAVO TOPMed germline variants
filtered_vcf_bravo_pass="${filtered_vcf_filter_pass//.vcf.gz/_bravo_pass.vcf.gz}"
if [[ ! -f "${filtered_vcf_bravo_pass}" ]]; then
    bcftools view "${filtered_vcf_filter_pass}" --targets-file "^${reference_directory}/$(basename ${bravo_variants})" --output "${filtered_vcf_bravo_pass}"
    tabix --force --preset vcf "${filtered_vcf_bravo_pass}"
fi

rsync -Prhltv ${filtered_vcf_bravo_pass} ${final_output_directory}
rsync -Prhltv ${filtered_vcf_bravo_pass}.tbi ${final_output_directory}
