# snakemake -s snakefile --profile folder/ --use-envmodules
# ask for 1Gb, 1CPU, 48h
import os

configfile: "/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/scripts/config.yaml"
#bams=["100_19965_20", "101_19226_20", "102_19843_20", "10_30_20","103_19902_20"]

chr_reg={}
chrs=[22]

def numberOfRegions():
    for chr in chrs:
        reg_file_dir=("chr_pos_test/chr%i" % chr) #change if regions file directory is changed.
        files = [f for f in os.listdir(reg_file_dir) if os.path.isfile(os.path.join(reg_file_dir, f))]
        nr_files=len(files)
        chr_reg[chr]=nr_files

'''
def output_files():
    numberOfRegions()
    out_list=[]
    for chr in chrs:
        for i in range(1, chr_reg[chr]):
            out_dir=("apply_var_recal/splits/chr%i/chr%i_reg%i.recalibrated.vcf.gz" % (chr,chr,i))
            out_list.append(out_dir)
    print(" ".join(out_list))
    return(" ".join(out_list))
'''

SAMPLES, = glob_wildcards("poles_example/{sample}.bam")
chromosomes, = glob_wildcards("chr_pos_test/chr{chr_id}/")


numberOfRegions() #number of splits per chromosome
out_files=[]
for chr in chrs:
    _out_per_chr=expand("apply_var_recal/splits/chr{c}/chr{c}_reg{i}.recalibrated.vcf.gz", c=chr, i=range(1,chr_reg[chr]))
    out_files.extend(_out_per_chr)

rule all:
    input:
        #lambda wildcards: output_files()
        out_files


#Call germline SNPs and indels via local re-assembly of haplotypes.
#capable of calling SNPs via local de-novo assembly of haplotypes in an active region.
#whenever the program encounters a region showing signs of variation, 
#it discards the existing mapping information and completely reassembles the reads in that region
rule HaplotypeCaller:
    input:
        #"poles_example/{SAMPLES}.bam"
        #config["samples"]
        "poles_example/{SAMPLES}.bam"
    output:
        "gvcf/chr{c}/{SAMPLES}.g.vcf.gz"
    log:
        "logs/HaplotypeCaller/chr{c}/{SAMPLES}.log"
    params:
        ref=config["ref"],
        chr_pos=config["chr_pos"]
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
        -L {params.chr_pos} \
        -ERC GVCF \
        -I {input} \
        --output-mode EMIT_ALL_ACTIVE_SITES \
        --native-pair-hmm-threads {resources.threads} \
        --min-pruning 2 --force-call-filtered-alleles true --alleles {params.chr_pos}
        """
'''
rule IndexVcf:
    input:
        "gvcf/chr{c}/{SAMPLES}.g.vcf.gz"
    output:
        "gvcf/chr{c}/{SAMPLES}.g.vcf.gz.tbi"
    envmodules:
        "bcftools"
    resources:
        mem='1g',
        time='00:30:0',
        threads=1
    shell:
        r"""
            bcftools index -t -f {input} -o {output}
        """
'''

rule GVCFSplit:
    input:
        gvcf="gvcf/chr{c}/{SAMPLES}.g.vcf.gz",
        reg="chr_pos_test/chr{c}/chr{c}_reg{i}.txt"
    output:
        "gvcf/splits/chr{c}/{SAMPLES}_reg{i}.g.vcf.gz"
    log:
        "gvcf/GVCFSplit/chr{c}/{SAMPLES}_reg{i}.log"
    envmodules:
        "bcftools"
    resources:
        mem='1g',
        time='00:30:0',
        threads=1
    shell:
        r"""
            bcftools +split {input.gvcf} -Oz -R {input.reg} > {output}
        """
'''
list_creater_input={}
for chr in chrs:
    list_creater=[]
    _list_per_chr=expand("gvcf/splits/chr{c}/{sample}_reg{i}.g.vcf.gz", c=chr, i=range(1,chr_reg[chr]), sample=SAMPLES)
    list_creater.extend(_list_per_chr)
    list_creater_input[chr]=list_creater
'''
def list_creater(wildcards):
    list_creater=[]
    _list_per_chr=expand("gvcf/splits/chr{c}/{sample}_reg{i}.g.vcf.gz", c=wildcards, i=range(1,chr_reg[chr]), sample=SAMPLES)
    list_creater.extend(_list_per_chr)
    return(list_creater)

#Creates a list of g.vcf.gz files created in the previous step, to be used in CombineGCVFs rule.
rule ListCreater:
    input:
        #"gvcf/splits/chr{c}/{SAMPLES}_reg{i}.g.vcf.gz"
        #("gvcf/splits/chr{c}/{SAMPLES}_reg{i}.g.vcf.gz")
        lambda wildcards: list_creater({wildcards.c})
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
        "gvcf/splits/chr{c}/combined_reg{i}.g.vcf.gz"
    params:
        ref=config["ref"],
        #int_1000=config["pos_1000"]
    envmodules:
        "gatk"
    resources:
        mem='20g',
        time='24:0:0',
        threads=8
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
        "gvcf/splits/chr{c}/combined_reg{i}.g.vcf.gz"
    output:
        "gvcf/splits/chr{c}/chr{c}_reg{i}.vcf.gz"
    params:
        ref=config["ref"],
        #int_1000=config["pos_1000"]
        chr_pos=config["chr_pos"]
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
        -V {input}  -O {output} \
        -L {params.chr_pos} -all-sites
        """ 

#Build a recalibration model to score variant quality for filtering purposes
#performs the first pass in Variant Quality Score Recalibration (VQSR). Builds the model to be used in ApplyVQSR.
#Developed adaptively based on "true sites" provided as input, then adaptive error model can be applied to known
#and novel variations to figure out the probability of each call to be true.
rule VariantRecalibrator:
    input:
        gvcf="gvcf/splits/chr{c}/chr{c}_reg{i}.vcf.gz",
        hap_map="/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/References/hapmap_3.3.hg38.vcf.gz",
        omni="/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/References/1000G_omni2.5.hg38.vcf.gz",
        TG="/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/References/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="/gpfs/space/GI/references/annotation_references/dbSNP/155/GCF_000001405.39_fixed_chr_names.gz"
    output:
        recal="var_recal/splits/chr{c}/chr{c}_reg{i}_sitesonly_AP.recal",
        tranch="var_recal/splits/chr{c}/chr{c}_reg{i}_sitesonly_AP.tranches",
        r_script="var_recal/splits/chr{c}/chr{c}_reg{i}_sitesonly_AP.plots.R"
    params:
        ref=config["ref"],
        #int_1000=config["pos_1000"],
        chr_pos=config["chr_pos"]
    envmodules:
        "gatk",
        "any/R/4.0.3"
    resources:
        mem='2g',
        time='24:0:0',
        threads=1       
    log:
        "logs/VariantRecalibrator/splits/chr{c}/chr{c}_reg{i}_sitesonly.log"
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
        -L {params.chr_pos} \
        --trust-all-polymorphic
        """   


#Second phase of VQSR, apply a score cutoff to filter variants based on a recalibration table. 
rule ApplyVQSR:
    input:
        gvcf="gvcf/splits/chr{c}/chr{c}_reg{i}.vcf.gz",
        recal="var_recal/splits/chr{c}/chr{c}_reg{i}_sitesonly_AP.recal",
        tranch="var_recal/splits/chr{c}/chr{c}_reg{i}_sitesonly_AP.tranches"
    log:
         "logs/ApplyVQSR/splits/chr{c}_reg{i}.recalibrated.log"
    envmodules:
        "gatk"
    resources:
        mem='5g',
        time='24:0:0',
        threads=4
    output:
        "apply_var_recal/splits/chr{c}/chr{c}_reg{i}.recalibrated.vcf.gz"
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
