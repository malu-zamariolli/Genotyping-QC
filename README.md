# Quality Control (QC)

# Usage

This Github contains the workflow pipeline that was used to perform quality control on genotype data. The workflow is organized as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. 

To run the pipeline with snakemake:
```console
cd  QC/
mamba activate snakemake
snakemake --use-conda --cores 1 -s snakefile
# number of cores can be modified accordingly
```

# Input data
For this study, we need:
- Genotype Data: plink binary files.

## Genotype Data 
Plink binary files (.bed, .bim, .fam) should be named in order to be recognized by defined wildcards in snakemake:

*{dataset}.bed -> EUR.bed*

*{dataset}.bim -> EUR.bim*

*{dataset}.fam -> EUR.fam*

# Config file

Threshold values for quality control (maf, geno, mind, hwe, pihat); genome-build and number of PCAs to be calculated can be defined in the *config.yaml* file 

# Tools

- [Plink v1.9](https://www.cog-genomics.org/plink/) can be used from environment created in the snakemake pipeline
    *../envs/plink.yaml*

# References
Some steps in this pipeline were adapted from [PRS-tutorial](https://choishingwan.github.io/PRS-Tutorial/prsice/).
