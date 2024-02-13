# snakemake -s snakefile --profile folder/ --use-envmodules
# ask for 1Gb, 1CPU, 48h
import os
from glob import glob

configfile: "/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/scripts/config.yaml"
#bams=["100_19965_20", "101_19226_20", "102_19843_20", "10_30_20","103_19902_20"]

chr_reg={}
chrs=list(range(1,23))

def numberOfRegions():
    for chr in chrs:
        reg_file_dir=("chr_pos_test/chr%i" % chr) #change if regions file directory is changed.
        files = [f for f in os.listdir(reg_file_dir) if os.path.isfile(os.path.join(reg_file_dir, f))]
        nr_files=len(files)
        chr_reg[chr]=nr_files



SAMPLES, = glob_wildcards("poles_example_full/{sample}.bam")
#chromosomes, = glob_wildcards("chr_pos_test/chr{chr_id}/")


numberOfRegions() #number of splits per chromosome
'''
out_files=[]
#expanding the output files before adding them into the rule all.
for chr in chrs:
    _out_per_chr=expand("apply_var_recal/splits/HRC_EST_POL.recalibrated.vcf.gz", c=chr, i=range(1,chr_reg[chr]+1))
    out_files.extend(_out_per_chr)
'''
'''
rule all:
    input:
        #lambda wildcards: output_files()
        #expand("apply_var_recal/chr{c}/HRC_EST_POL.recalibrated.vcf.gz", c=chrs)
        "norm_var_recal/HRC_EST_POL.recalibrated.encoded.norm.vcf.gz"
'''

rule all:
    input:
        "norm_var_recal/HRC_EST_POL.recalibrated.encoded.norm.snps.vcf.gz"


#Call germline SNPs and indels via local re-assembly of haplotypes.
#capable of calling SNPs via local de-novo assembly of haplotypes in an active region.
#whenever the program encounters a region showing signs of variation, 
#it discards the existing mapping information and completely reassembles the reads in that region

rule HaplotypeCaller:
    input:
        #"poles_example/{SAMPLES}.bam"
        #config["samples"]
        bam="poles_example/{SAMPLES}.bam",
        region="/gpfs/space/home/alacamli/HRC_EST_POL/Ref_Enrichment/HRC_chr_pos/chr{c}_autosomal_filtered_reheader.vcf.gz"
    output:
        "gvcf/chr{c}/{SAMPLES}.g.vcf.gz"
    log:
        "logs/HaplotypeCaller/chr{c}/{SAMPLES}.log"
    params:
        ref=config["ref"]
        #chr_pos=config["chr_pos"]
    benchmark:
        "benchmarks/HaplotypeCaller/chr{c}_{SAMPLES}.benchmark.txt"
    envmodules:
        "gatk"
    resources:
        mem='5g',
        time='24:0:0',
        threads=1
    shell:
        r"""
        gatk --java-options -Xmx{resources.mem} HaplotypeCaller --add-output-vcf-command-line \
        -R {params.ref} \
        -O {output} \
        -L {input.region} \
        -ERC GVCF \
        -I {input.bam} \
        --output-mode EMIT_ALL_ACTIVE_SITES \
        --create-output-variant-index true \
        --native-pair-hmm-threads {resources.threads} \
        --min-pruning 2 --force-call-filtered-alleles true --alleles {input.region}
        """


rule GVCFSplit:
    input:
        gvcf="gvcf/chr{c}/{SAMPLES}.g.vcf.gz",
        reg="chr_pos_test/chr{c}/chr{c}_reg{i}.txt"
    output:
        "gvcf/splits/chr{c}/{SAMPLES}_reg{i}.g.vcf.gz"
    log:
        "logs/GVCFSplit/chr{c}/{SAMPLES}_reg{i}.log"
    benchmark:
        "benchmarks/GVCFSplit/chr{c}/chr{c}_{SAMPLES}_reg{i}.benchmark.txt"
    envmodules:
        "bcftools"
    resources:
        mem='1g',
        time='00:30:0',
        threads=1
    shell:
        r"""
            bcftools view {input.gvcf} -Oz -R {input.reg} -o {output}
            tabix -p vcf {output}
        """

def list_creater(chr, reg):
    list_creater=[]
    _list_per_chr=expand("gvcf/splits/chr{c}/{sample}_reg{i}.g.vcf.gz", c=chr, i=reg, sample=SAMPLES)
    list_creater.extend(_list_per_chr)
    return(list_creater)

#Creates a list of g.vcf.gz files created in the previous step, to be used in CombineGCVFs rule.
rule ListCreater:
    input:
        #"gvcf/splits/chr{c}/{SAMPLES}_reg{i}.g.vcf.gz"
        #("gvcf/splits/chr{c}/{SAMPLES}_reg{i}.g.vcf.gz")
        lambda wildcards: list_creater({wildcards.c}, {wildcards.i})
    output:
        "gvcf/splits/chr{c}/chr{c}_reg{i}_gvcf.list"
    shell:
        r"""
        ls {input} > {output}
        """


#Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file.
rule CombineGVCFs:
    input:
        #expand("gvcf/chr{c}/{sample}.g.vcf.gz", sample=bams, c=chrs)
        "gvcf/splits/chr{c}/chr{c}_reg{i}_gvcf.list"
    output:
        "gvcf/splits/chr{c}/combined/combined_reg{i}.g.vcf.gz"
    params:
        ref=config["ref"],
        #int_1000=config["pos_1000"]
    benchmark:
        "benchmarks/CombineGVCFs/chr{c}/chr{c}_reg{i}.benchmark.txt"
    envmodules:
        "gatk"
    resources:
        mem='20g',
        time='24:0:0',
        threads=1
    log:
        "logs/CombineGVCFs/splits/chr{c}_reg{i}.log"
    shell:
        r"""
        gatk --java-options -Xmx{resources.mem} CombineGVCFs \
        -R {params.ref} \
        -O {output} -V {input}
        """

#Perform joint genotyping on one or more samples pre-called with HaplotypeCaller.

rule GenotypeGVCFs:
    input:
        combined="gvcf/splits/chr{c}/combined/combined_reg{i}.g.vcf.gz",
        region="/gpfs/space/home/alacamli/HRC_EST_POL/Ref_Enrichment/HRC_chr_pos/chr{c}_autosomal_filtered_reheader.vcf.gz"
    output:
        "GenotypeGVCFs/splits/chr{c}/chr{c}_reg{i}.vcf.gz"
    params:
        ref=config["ref"],
        #int_1000=config["pos_1000"]
        chr_pos=config["chr_pos"]
    benchmark:
        "benchmarks/GenotypeGVCFs/chr{c}/chr{c}_reg{i}.benchmark.txt"
    envmodules:
        "gatk"
    resources:
        mem='20g',
        time='24:0:0',
        threads=1
    log:
        "logs/GenotypeGVCFs/splits/chr{c}_reg{i}.log"
    shell:
        r"""
        gatk --java-options -Xmx{resources.mem} GenotypeGVCFs \
        -R {params.ref} \
        -V {input.combined}  -O {output} \
        -L {input.region} -all-sites
        """

concat_input_files=[]
#expanding the output files before adding them into the rule all.
for ch in chrs:
    _in_per_chr=expand("GenotypeGVCFs/splits/chr{c}/chr{c}_reg{i}.vcf.gz", c=ch, i=range(1,chr_reg[ch]+1))
    concat_input_files.extend(_in_per_chr)

rule bcftoolsConcat:
    input:
        concat_input_files
    output:
        "bcftoolsConcat/HRC_EST_POL.vcf.gz"
    benchmark:
        "benchmarks/bcftoolsConcat/HRC_EST_POL.txt" #CHANGE WHEN RUNNING FOR WHOLE GENOME
    resources:
        mem='20g',
        time='24:0:0',
        threads=1        
    log:
        "logs/bcftoolsConcat/HRC_EST_POL.log"
    envmodules:
        "bcftools"
    shell:
        r"""
            bcftools concat -Oz -o {output} {input}
            tabix -p vcf {output}
        """

#Build a recalibration model to score variant quality for filtering purposes
#performs the first pass in Variant Quality Score Recalibration (VQSR). Builds the model to be used in ApplyVQSR.
#Developed adaptively based on "true sites" provided as input, then adaptive error model can be applied to known
#and novel variations to figure out the probability of each call to be true.
rule VariantRecalibrator:
    input:
        gvcf="bcftoolsConcat/HRC_EST_POL.vcf.gz",
        hap_map="/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/References/hapmap_3.3.hg38.vcf.gz",
        omni="/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/References/1000G_omni2.5.hg38.vcf.gz",
        TG="/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/References/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="/gpfs/space/home/alacamli/HRC_EST_POL/Ref_Enrichment/References/Homo_sapiens_assembly38.dbsnp138.vcf7"
    output:
        recal="var_recal/HRC_EST_POL_sitesonly_AP.recal",
        tranch="var_recal/HRC_EST_POL_sitesonly_AP.tranches",
        r_script="var_recal/HRC_EST_POL_sitesonly_AP.plots.R"
    params:
        ref=config["ref"],
        #int_1000=config["pos_1000"],
        chr_pos=config["chr_pos"]
    envmodules:
        "gatk",
        "any/R/4.0.3"
    benchmark:
        "benchmarks/VariantRecalibrator/HRC_EST_POL.txt"
    resources:
        mem='30g',
        time='24:0:0',
        threads=1       
    log:
        "logs/VariantRecalibrator/HRC_EST_POL_sitesonly.log"
    shell:
        r"""
        gatk VariantRecalibrator \
        -R {params.ref} \
        -V {input.gvcf} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hap_map} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.TG} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O {output.recal} \
        --tranches-file {output.tranch} \
        --rscript-file {output.r_script} \
        --dont-run-rscript \
        -L {params.chr_pos} \
        --trust-all-polymorphic
        """   


#Second phase of VQSR, apply a score cutoff to filter variants based on a recalibration table. 
rule ApplyVQSR:
    input:
        gvcf="bcftoolsConcat/HRC_EST_POL.vcf.gz",
        recal="var_recal/HRC_EST_POL_sitesonly_AP.recal",
        tranch="var_recal/HRC_EST_POL_sitesonly_AP.tranches"
    log:
         "logs/ApplyVQSRHRC_EST_POL.recalibrated.log"
    envmodules:
        "gatk"
    benchmark:
        "benchmarks/ApplyVQSRHRC_EST_POL.txt"    
    resources:
        mem='10g',
        time='24:0:0',
        threads=1
    output:
        "apply_var_recal/HRC_EST_POL.recalibrated.vcf.gz"
    shell:
        r"""
        gatk --java-options -Xmx{resources.mem} ApplyVQSR \
        -V {input.gvcf} \
        --recal-file {input.recal} \
        --tranches-file {input.tranch} \
        --truth-sensitivity-filter-level 99.7 \
        --create-output-variant-index true \
        -mode SNP -O {output}
        """
    
rule encodeMissingness:
    input:
        "apply_var_recal/HRC_EST_POL.recalibrated.vcf.gz"
    output:
        "encoded_var_recal/HRC_EST_POL.recalibrated.encoded.vcf.gz"
    envmodules:
        "bcftools"
    log:
        "logs/encodeMissingness/HRC_EST_POL.recalibrated.encoded.log"
    benchmark:
        "benchmarks/encodeMissingness/HRC_EST_POL.recalibrated.txt"
    resources:
        mem='1g',
        time='24:0:0',
        threads=1
    shell:
        r"""
        bcftools +fill-tags {input} -Oz -o {output} -- -t F_MISSING
        tabix -p vcf {output}
        """

rule NormVCF:
    input:
        "encoded_var_recal/HRC_EST_POL.recalibrated.encoded.vcf.gz"
    output:
        "norm_var_recal/HRC_EST_POL.recalibrated.encoded.norm.all.vcf.gz"
    envmodules:
        "bcftools"
    log:
        "logs/normVCF/HRC_EST_POL.recalibrated.encoded.norm.all.log"
    benchmark:
        "benchmarks/normVCF/HRC_EST_POL.recalibrated.encoded.norm.all.txt"
    resources:
        mem='5g',
        time='24:0:0',
        threads=1
    shell:
        r"""
        bcftools norm {input} -Oz -o {output} -m -
        tabix -p vcf {output}
        """
rule filterSNPS:
    input:
        "norm_var_recal/HRC_EST_POL.recalibrated.encoded.norm.all.vcf.gz"
    output:
        "norm_var_recal/HRC_EST_POL.recalibrated.encoded.norm.snps.vcf.gz"
    envmodules:
        "bcftools"
    log:
        "logs/filterSNPs/HRC_EST_POL.recalibrated.encoded.norm.snps.log"
    benchmark:
        "benchmarks/filterSNPs/HRC_EST_POL.recalibrated.encoded.norm.snps.txt"
    resources:
        mem='5g',
        time='24:0:0',
        threads=1
    shell:
        r"""
        bcftools norm {input} -Oz -o {output} -m -
        tabix -p vcf {output}
        """


