
rule all:
	input:
		expand("SNPFiltered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.snp_filtered.stats.txt", f=["PASS","HRC_sites"]),
		expand("Ind_Filtered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.stats.txt", f=["PASS","HRC_sites"]),
		expand("Ind_stats/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.stats.txt", f=["PASS","HRC_sites"])

rule HRCSites:
	input:
		ref="norm_var_recal/HRC_EST_POL.recalibrated.encoded.norm.snps.vcf.gz",
		hrc="/gpfs/space/GI/ebc_data/projects/HRC_EST_POL/Ref_Enrichment/HRC_chr_pos/HRC_WG_autosomal_filtered_reheader.sort.vcf.gz"
	output:
		directory("HRC_sites/")
	envmodules:
		"bcftools"
	benchmark:
		"benchmarks/hrc_sites/HRC_EST_POL.recalibrated.encoded.norm.snps.HRC_sites.txt"
	log:
		"logs/hrc_sites/HRC_EST_POL.recalibrated.encoded.norm.snps.HRC_sites.log"
	resources:
		mem='5g',
		time='24:0:0',
		threads=1
	shell:
		r"""
			mkdir {output}
			bcftools isec -c none -n=2 -w 2 {input.hrc} {input.ref} -p {output}
			mv {output}/0001.vcf {output}/HRC_EST_POL.recalibrated.encoded.norm.snps.HRC_sites.vcf
			bgzip -c {output}/HRC_EST_POL.recalibrated.encoded.norm.snps.HRC_sites.vcf > {output}/HRC_EST_POL.recalibrated.encoded.norm.snps.HRC_sites.vcf.gz
			tabix -p vcf {output}/HRC_EST_POL.recalibrated.encoded.norm.snps.HRC_sites.vcf.gz
		"""
#


rule FilterPASS:
	input:
		"HRC_sites/HRC_EST_POL.recalibrated.encoded.norm.snps.HRC_sites.vcf.gz"
	output:
		"HRC_sites/HRC_EST_POL.recalibrated.encoded.norm.snps.PASS.vcf.gz"
	envmodules:
		"bcftools"
	log:
		"logs/FilterPASS/HRC_EST_POL.recalibrated.encoded.norm.snps.PASS.log"
	benchmark:
		"benchmarks/FilterPASS/HRC_EST_POL.recalibrated.encoded.norm.snps.PASS.txt"
	resources:
		mem='1g',
		time='24:0:0',
		threads=1
	shell:
		r"""
		bcftools filter -i 'FILTER="PASS"' {input} -Oz -o {output}
		tabix -p vcf {output}
		"""
rule IndStats:
	input:
		"HRC_sites/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.vcf.gz"
	output:
		"Ind_stats/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.stats.txt"
	envmodules:
		"bcftools"
	log:
		"logs/IndStats/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.log"
	benchmark:
		"benchmarks/IndStats/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.stats.txt"
	resources:
		mem='1g',
		time='24:0:0',
		threads=1
	shell:
		r"""
		bcftools stats -s- {input} > {output}
		"""

rule Ind_Missingness:
	input:
		"HRC_sites/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.vcf.gz"
	output:
		samples="Ind_stats/HRC_EST_POL.samples.{f}.txt",
		missingness="Ind_stats/HRC_EST_POL.missingness.{f}.txt",
		ind_missingness="Ind_stats/HRC_EST_POL.ind_missingness.{f}.txt",
		filtering_samples="Ind_stats/HRC_EST_POL.filtering_samples.{f}.txt"
	envmodules:
		"bcftools"
	log:
		"logs/Ind_Missingness/HRC_EST_POL.ind_missingness.{f}.log"
	benchmark:
		"benchmarks/Ind_Missingness/HRC_EST_POL.ind_missingness.{f}.txt"
	resources:
		mem='1g',
		time='3:0:0',
		threads=1
	shell:
		r"""
			./ind_stats_missingness.sh {input} {output.samples} {output.missingness} {output.ind_missingness} {output.filtering_samples}
		"""

rule Ind_Filtering:
	input:
		vcf_file="HRC_sites/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.vcf.gz",
		filtering_samples="Ind_stats/HRC_EST_POL.filtering_samples.{f}.txt"
	output:
		"Ind_Filtered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.vcf.gz"
	envmodules:
		"bcftools"
	log:
		"logs/Ind_Filtering/HRC_EST_POL.ind_filtered.{f}.log"
	benchmark:
		"benchmarks/Ind_Filtering/HRC_EST_POL.ind_filtered.{f}.txt"
	resources:
		mem='1g',
		time='3:0:0',
		threads=1
	shell:
		r"""
			bcftools view -S ^{input.filtering_samples} -Oz -o {output} {input.vcf_file}
		"""

rule StatsAfterIndFilter:
	input:
		"Ind_Filtered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.vcf.gz"
	output:
		"Ind_Filtered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.stats.txt"
	envmodules:
		"bcftools"
	log:
		"logs/StatsAfterIndFilter/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.stats.log"
	benchmark:
		"benchmarks/StatsAfterIndFilter/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.stats.txt"
	resources:
		mem='1g',
		time='3:0:0',
		threads=1
	shell:
		r"""
			bcftools stats {input} > {output}
		"""

rule SnpFiltering:
	input:
		"Ind_Filtered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.vcf.gz"
	output:
		"SNPFiltered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.snp_filtered.vcf.gz"
	envmodules:
		"bcftools"
	log:
		"logs/SnpFiltering/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.snp_filtered.log"
	benchmark:
		"benchmarks/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.snp_filtered.txt"
	resources:
		mem='1g',
		time='3:0:0',
		threads=1
	shell:
		r"""
			bcftools view -i 'F_MISSING > 0.05' -Oz -o {output} {input}
		"""

rule FinalStats:
	input:
		"SNPFiltered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.snp_filtered.vcf.gz"
	output:
		"SNPFiltered/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.snp_filtered.stats.txt"
	envmodules:
		"bcftools"
	log:
		"logs/FinalStats/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.snp_filtered.stats.log"
	benchmark:
		"benchmarks/HRC_EST_POL.recalibrated.encoded.norm.snps.{f}.ind_filtered.snp_filtered.stats.benchmark.txt"
	resources:
		mem='1g',
		time='3:0:0',
		threads=1
	shell:
		r"""
			bcftools stats {input} > {output}
		"""