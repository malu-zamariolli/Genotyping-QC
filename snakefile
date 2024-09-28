#----------------------- Genotyping Quality Control --------------------------------#
# A. Config file
configfile:"config.yaml"

# B. Define wildcard
DATASETS,FORMATS = glob_wildcards("../Binary/{dataset}.{ext}")

# C. Rule all
rule all:
     input:
          expand(["QC_reports/{dataset}_QC_report.pdf"],dataset=DATASETS),
          expand(["QC_reports/{dataset}_sex_check.pdf"], dataset=DATASETS),
          expand(["Binary_final/{dataset}_QCed_final.{ext}"],  dataset=DATASETS, ext=FORMATS),
          expand(["PCA/{dataset}_plotPCA.pdf"], dataset=DATASETS),
          expand(["Binary_final/{dataset}_N_SNPs_variants_final.txt"], dataset=DATASETS)
          # Rule all only listing files that are output and not input to another rule


# 1. Remove indels and SNPs with more than 2 alleles
rule biallelic:
     input:
          multiext("../Binary/{dataset}", ".bed", ".bim", ".fam")
     output:
          temp(multiext("QC_binary_temp/{dataset}_biallelic", ".bed", ".bim", ".fam"))
     params:
        a1 = "../Binary/{dataset}",
        a2 = "QC_binary_temp/{dataset}_biallelic"
     conda: "./envs/plink.yaml"
     shell: 
          """
          plink --bfile {params.a1} --snps-only just-acgt --make-bed --out {params.a2}
          """
# 2. QC reports
## 2. A) Check MAF, missing SNPs, Hardy-Weinberg Equilibrium, heterozigosity
rule QC_report:
     input:
          multiext("QC_binary_temp/{dataset}_biallelic", ".bed", ".bim", ".fam")
     output:
          multiext("QC_reports/{dataset}_QC_reports", ".frq", ".lmiss", ".imiss", ".hwe", ".het")
     params:
          a2 = "QC_binary_temp/{dataset}_biallelic",
          a3 = "QC_reports/{dataset}_QC_reports"
     conda: "./envs/plink.yaml"
     shell:
          """
          plink --bfile {params.a2} --freq --missing --hardy --het --out {params.a3} 
          """

## 2. B) Plot MAF, Individual missingness, Locus missingness, Hardy_Weinberg p-value
rule QCreport_plot:
     input:
          imiss = "QC_reports/{dataset}_QC_reports.imiss",
          lmiss = "QC_reports/{dataset}_QC_reports.lmiss", 
          maf = "QC_reports/{dataset}_QC_reports.frq",
          het = "QC_reports/{dataset}_QC_reports.het",
          hwe = "QC_reports/{dataset}_QC_reports.hwe"
     output:
          plots = "QC_reports/{dataset}_QC_report.pdf" 
     script:
          "scripts/qc_reports_plot.R"
  
## 3. QC: Filter SNP missingness (geno)
rule geno:
     input:
          multiext("QC_binary_temp/{dataset}_biallelic", ".bed", ".bim", ".fam")
     output:
          temp(multiext("QC_binary_temp/{dataset}_QCedSNP", ".bed", ".bim", ".fam"))
     params:
          a2 = "QC_binary_temp/{dataset}_biallelic",
          a4 = "QC_binary_temp/{dataset}_QCedSNP",
          geno = config["geno"]
     conda: "./envs/plink.yaml"
     shell:
          """
          plink --bfile {params.a2} --geno {params.geno} --make-bed --out {params.a4}
          """
# --geno: # filters out all variants with missing call rates exceeding the provided value

## 3. QC: Filter individual missingness (mind) and MAF, Hardy-Weinberg
rule QC:
     input:
          multiext("QC_binary_temp/{dataset}_QCedSNP", ".bed", ".bim", ".fam")
     output:
          temp(multiext("QC_binary_temp/{dataset}_QCed", ".bed", ".bim", ".fam"))
     params:
          a2 = "QC_binary_temp/{dataset}_QCedSNP",
          a4 = "QC_binary_temp/{dataset}_QCed",
          maf = config["maf"],
          mind = config["mind"],
          hwe = config["hwe"]
     conda: "./envs/plink.yaml"
     shell:
          """
          plink --bfile {params.a2} --mind {params.mind} --maf {params.maf} --hwe {params.hwe} --make-bed --out {params.a4}
 
          """
# --mind is being performed after geno to protect for loosing too many samples in the cohort
# --hwe filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold
     # Be aware: with case/control data, cases and missing phenotypes are normally ignored (--hwe)

## 4. QC: Heterozigosity
# Heterozigsity tests should be done on a set of SNPs that are not highly correlated
rule prune:
    input:
          multiext("QC_binary_temp/{dataset}_QCed", ".bed", ".bim", ".fam")
    output:
          multiext("lists/{dataset}_QCed", ".prune.in", ".prune.out")
    params:
        a4 = "QC_binary_temp/{dataset}_QCed",
        a5 = "lists/{dataset}_QCed"
    conda: "./envs/plink.yaml"
    shell:
        """
        plink \
        --bfile {params.a4} \
        --indep-pairwise 200 50 0.25 \
        --out {params.a5}
        """
rule heterozygosity:
     input:
          multiext("QC_binary_temp/{dataset}_QCed", ".bed", ".bim", ".fam"),
          prune_in = "lists/{dataset}_QCed.prune.in"
     output: 
          "QC_reports/{dataset}_QC_reports_pruned.het"
     params:
        a4 = "QC_binary_temp/{dataset}_QCed",
        a6 = "QC_reports/{dataset}_QC_reports_pruned"
     conda: "./envs/plink.yaml"
     shell:
        """
        plink \
        --bfile {params.a4} \
        --extract {input.prune_in} \
        --het \
        --out {params.a6}
        """

rule filter_het:
     input:
        het_qc = "QC_reports/{dataset}_QC_reports_pruned.het"
     output:
        het_out = "lists/{dataset}_het_valid_sample"
     script:
        "scripts/heterozygosity.R"

## 5. QC: Sex-check
rule sexcheck:
     input:
          multiext("QC_binary_temp/{dataset}_QCed", ".bed", ".bim", ".fam"),
          prune_in = "lists/{dataset}_QCed.prune.in",
          het_out = "lists/{dataset}_het_valid_sample"
     output:
          sexcheck = "lists/{dataset}_QCed.sexcheck"
     params:
          a9 = "QC_binary_temp/{dataset}_QCed",
          a7 = "lists/{dataset}_QCed"
     conda: "./envs/plink.yaml"        
     shell:
        """
       plink \
            --bfile {params.a9} \
            --extract {input.prune_in} \
            --keep {input.het_out} \
            --check-sex \
            --out {params.a7}
        """
# In this step we are only keeping the samples after heterozigozity filtering
# Females - F <0.2 | Males F > 0.8

rule filter_sex:
    input:
        valid = "lists/{dataset}_het_valid_sample",
        sexcheck = "lists/{dataset}_QCed.sexcheck"
    output:
        filtered_sex = "lists/{dataset}_het_sexcheck_valid_sample",
        plot = "QC_reports/{dataset}_sex_check.pdf"
    script:
        "scripts/filter_sexcheck.R"

## 6. Remove chromosomes 
rule remove_chr:
     input:
          multiext("QC_binary_temp/{dataset}_QCed", ".bed", ".bim", ".fam"),
          valid_sample = "lists/{dataset}_het_sexcheck_valid_sample"
     output:
          temp(multiext("QC_binary_temp/{dataset}_QCed_noXY", ".bed", ".bim", ".fam"))
     params:
          a9 = "QC_binary_temp/{dataset}_QCed",
          a10 = "QC_binary_temp/{dataset}_QCed_noXY"
     conda: "./envs/plink.yaml"
     shell:
          """
          plink  \
               --bfile {params.a9}  \
               --not-chr 0,23-26  \
               --keep {input.valid_sample} \
               --make-bed  \
               --keep-allele-order  \
               --out {params.a10}
          """
# (0: unknown, 23: X, 24: Y, 25: XY PAR, 26: mitochondrial)
# input binary from before splitsex rule (since that rule splitsex only keeps pruned SNPs in the output binary - only used for sexcheck)

## 7. Identity by descent
rule ibd:
     input:
          multiext("QC_binary_temp/{dataset}_QCed_noXY", ".bed", ".bim", ".fam"),
          prune_in = "lists/{dataset}_QCed.prune.in"
     output:
          "lists/{dataset}_ibd_calculation.genome"
     params:
          a10 = "QC_binary_temp/{dataset}_QCed_noXY",
          a11 = "lists/{dataset}_ibd_calculation",
          pihat = config["pihat"]
     conda: "./envs/plink.yaml"
     shell:
          """
          plink --bfile {params.a10} --extract {input.prune_in} --genome --min {params.pihat} --out {params.a11}
          """
# this calculation is not LD-aware, so it's good to do pruning before
# Pihat 0.125: Third-degree relatives (12.5% equal IBD) (first cousins)

rule relatedness_list:
     input:
          ibd = "lists/{dataset}_ibd_calculation.genome"
     output:
          rel_list = "lists/{dataset}_related_individuals_exclude.txt"
     shell:
          "awk '{{ print $3,$4 }}' {input.ibd} > {output.rel_list}"

rule remove_related:
     input:
          multiext("QC_binary_temp/{dataset}_QCed_noXY", ".bed", ".bim", ".fam"),
          rel_list = "lists/{dataset}_related_individuals_exclude.txt"
     output:
          multiext("Binary_final/{dataset}_QCed_final", ".bed", ".bim", ".fam"),
     params:
          a10 = "QC_binary_temp/{dataset}_QCed_noXY",
          a12 = "Binary_final/{dataset}_QCed_final"
     conda: "./envs/plink.yaml"
     shell:
          """
          plink --bfile {params.a10} --remove {input.rel_list} --make-bed --out {params.a12}
          """
## 8. PCA
rule pca:
     input:
          multiext("Binary_final/{dataset}_QCed_final", ".bed", ".bim", ".fam")
     output:
          eigenvec = "PCA/{dataset}_pca.eigenvec",
          eigenval = "PCA/{dataset}_pca.eigenval"
     params:
          a12 = "Binary_final/{dataset}_QCed_final",
          a13 = "PCA/{dataset}_pca",
          pca = config["pca"]
     conda: "./envs/plink.yaml"
     shell:
          "plink --bfile {params.a12} --pca {params.pca} header --out {params.a13}"

rule plot_pca:
     input:
          eigenvec = "PCA/{dataset}_pca.eigenvec"
     output:
          plot_pca = "PCA/{dataset}_plotPCA.pdf"
     script:
           "scripts/plot_PCA.R"

## 9. Summary report
rule final_report:
     input:
          initial_bim = "../Binary/{dataset}.bim",
          initial_fam = "../Binary/{dataset}.fam",
          final_bim = "Binary_final/{dataset}_QCed_final.bim",
          final_fam = "Binary_final/{dataset}_QCed_final.fam"
     output:
          report = "Binary_final/{dataset}_N_SNPs_variants_final.txt"
     script:
          "scripts/variants_individuals.R"