# Ref Enrichment

## Script Location
The snakemake script is located in scripts/Ref_Enrichment.smk, the bash script to submit the job is located in scripts/Ref_Enrichment.sh. The config.yaml file to pass the variables to the snakemake code is located in scritps/config.yaml
and the HPC configuration file is located in scripts/config_dir.

## Workflow

### HaplotypeCaller
  The workflow starts with GATK HaplotypeCaller function. The rule takes .bam files as input and is used to call potential variant sites per sample and save results in GVCF format.  
  ![image](https://github.com/SonofLastPanda/Ref_Enrichment/assets/41624456/56f2d0d8-f343-4fc1-968a-41c1d6788d7f)  

1. **Define Active Regions:** HaplotypeCaller determines which regions of the genome it needs to operate on, given by -R option.
2. **Determine haplotypes by re-assembly of the active region:** For each active region, the program builds a De Bruijn-like graph to identify the possible haplotypes in the data.
Then, the haplotypes are mapped against the reference haplotype to identify the potential variants.
3. **Determine likelihoods of the haplotypes given the read data:** For each active region, do pairwise alignment of each read to each haplotype with PairHMM algorithm, which produces a matrix of likelihoods of haplotypes given the read data, then the likelihoods
are marjinalized to obtain the likelihoods of alleles per read for each potentially variant site.
4. **Assign sample genotypes:**  Bayes' rule is applied for each potential variant, using the likelihoods of alleles given the read data to calculate the posterior likelihoods of each genotype per sample given the read data observed for that sample.
The most likely genotype is then assigned to the sample.

### GVCFSplit
  The GVCF files are splitted into regions of roughly 50k SNPs, using ``` bcftools view ``` with ``` -R ``` flag, to make the combination step to consume less memory. The splits are indexed afterwards with ```  bcftools index ```.

### ListCreater
  A rule to create the list of input files to be used in CombineGVCFs rule. The list is created for each region and chromosome usign all the samples.

### CombineGVCFs
  Combine the splits of per-sample GVCFs into multi-sample GVCF file.

### GenotypeGVCFs
  Perform joint genotyping on the files outputted by CombineGVCFs file

### VariantRecalibrator
  First pass of in Variant Quailty Score Recalibration (VQSR). Develops an adaptive error model based on "true sites", the error model then can be applied to known and novel variations to determine their probability of being true.

### ApplyVSQR
  Second phase VQSR. The model generated in the previous step is applied to genotype data, and each variant is marked as pass/fail in the output VCF data, according to the threshold given. 
The failed variants are not removed from the VCF unless specified, rather marked as filtered.

### Workflow Graph

![dag](https://github.com/SonofLastPanda/Ref_Enrichment/assets/41624456/635fcedf-1952-4b0f-b3fa-5fa97537de31)
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
 "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<!-- Generated by graphviz version 2.49.0 (0)
 -->
<!-- Title: snakemake_dag Pages: 1 -->


<!-- 7&#45;&gt;1 -->

</svg>



