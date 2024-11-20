This repository contains the pipeline used for calling passenger mutations from whole genome sequencing (WGS) data in the UK BioBank.  The dependencies are all available on bioconda and can be installed using [pixi](https://pixi.sh/latest/):
```
pixi global install -c bioconda -c conda-forge bcftools gatk4 htslib
```

The other components are standard Unix tools but you can optionally install newer versions with [pixi](https://pixi.sh/latest/) if desired:
```
pixi global install findutils grep parallel rsync sed util-linux
pixi global expose remove kill  # needed to avoid conflict with coreutils
pixi global install coreutils
```

It has 5 components (along with some old scripts) that are used in this order:

1. **mutect_wgs_scatter**: This step is run on the UK Biobank Research Analysis Platform. It runs the two steps in the pipeline with *GATK* that require access to the raw CRAM files: *Mutect2* and *GetPileupSummaries*. It "scatters" these operations across chromosomes so they can be run in parallel.
2. **mutect_wgs_gather**:  This step is run on our local SLURM-based cluster.  It aggregates the "scattered" output from *Mutect2* and *GetPileupSummaries* into single files per sample and then applies several annotation and filtering steps:
    - *LearnReadOrientationModel* (*GATK*) - learn a read orientation bias model to be used as an input for *FilterMutectCalls*
    - *CalculateContamination* (*GATK*) - calculate cross-sample contamination rate as an input for *FilterMutectCalls*
    - *FilterMutectCalls* (*GATK*) - annotation the `FILTER` column in the VCF with information about whether each variant call passes or fails various quality control filters
    - Remove variants which are not biallelic single nucleotide variants (SNVs) with *bcftools* - multiallelic SNVs and indels are far more likely to be artifacts than biallelic SNVs
    - Remove variants which are found in the TOPMed Bravo list of common germline variants - these variants are likely to be germline artifacts
3. **bcftools_merge**: This step is run on our local SLURM-based cluster.  It uses *bcftools* to merge VCFs from the previous step in groups of 100 samples.  This is more performant than trying to merge thousands of samples at once.  The merged VCFs are split by chromosome to facilitate parallelization in later steps.  The merged VCFs have also been stripped of all tags that are not needed for subsequent steps to reduce file size.
4. **bcftools_filter**: This step is run on our local SLURM-based cluster.  It uses *bcftools* to merge the VCFs from the previous step across into a single VCF per chromosome.  Allele counts are computed with the *+fill-tags* plugin from *bcftools*.  It then applies the following filters with `bcftools`:
    - Allelic depth of the alternate allele greater than or equal to 4: `AD[*:1]>=4`
    - Total depth less than or equal to 100: `DP<=100`
    - Allele frequency less than or equal to 0.4: `AF<=0.4`
    - Variant passes FilterMutectCalls filters: `FILTER='PASS'`
    - Variant is a singleton (allele count of 1): `INFO/AC=1`
    - Variant is still biallelic after merging across all samples
5. **bcftools_concat**: This step is run on our local SLURM-based cluster.  It uses *bcftools* to concatenate the filtered per chromosome VCFs from the previous step into a single VCF.  It then uses *bcftools* to count the number of singletons per sample for all remaining variants, as well as only C>T/T>C variants.  Unused code is also provided for removing variants in "difficult to sequence" intervals provided by Genome in a Bottle and removing variants in the mDust and superDups intervals marking repetitive and duplicated regions of the genome

Please note that for all steps that are run on our SLURM-based cluster, the "submit" scripts are primarily used for logic specific to that job schedule and would need modifications to work with other schedulers or cloud-based systems.  

Also note that all usage of `rsync` is only for transferring files to local high-performance storage and may not be relevant for all compute environments.  We do not recommend running the most disk-intensive steps on "archival" storage such as cloud-based AWS S3/Azure Blob Storage/Google Cloud Storage or local HPC systems such as LustreFS.  Instead, inputs should be copied to higher performance storage such as cloud-based AWS EBS/Azure Managed Disks/Google Cloud Persistent Disks or local HPC flash storage (NVMe preferred over SATA) should be used for these steps, after which outputs can be copied back to archival storage.

This workflow may eventually be ported to [nf-core](https://nf-co.re).  Nextflow was not used in the initial design due to the limited support for Nextflow on DNA Nexus, and the lack of support within Nextflow for SLURM job arrays, but will probably make the more complex workflow logic simpler to maintain.
