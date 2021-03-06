#!/usr/bin/python

import os
import sys
import yaml

configfile: "config.yaml"
print(config)

samples = [sample for sample in config['samples']]
print(samples)

if not os.path.isdir('analysis'):
    os.mkdir('analysis')

def get_outputs(config, samples):
    ls = list()
    ls.append('{0}.amb'.format(config['reference']['genome_fasta']))
    ls.append('{0}.ann'.format(config['reference']['genome_fasta']))
    ls.append('{0}.bwt'.format(config['reference']['genome_fasta']))
    ls.append('{0}.pac'.format(config['reference']['genome_fasta']))
    ls.append('{0}.sa'.format(config['reference']['genome_fasta']))
    ls.append('{0}.fai'.format(config['reference']['genome_fasta']))
    ls.append('{0}.dict'.format(config['reference']['genome_fasta'].rstrip('.fasta')))

    for sample in samples:
        #ls.append('analysis/fastqc/{0}/{0}_1_fastqc.zip'.format(sample))
        #ls.append('analysis/fastqc/{0}/{0}_2_fastqc.zip'.format(sample))

        ls.append('analysis/sickle/{0}/{0}_R1.fastq.gz'.format(sample))
        ls.append('analysis/sickle/{0}/{0}_R2.fastq.gz'.format(sample))
        ls.append('analysis/sickle/{0}/{0}_se.fastq.gz'.format(sample))

        #ls.append('analysis/statistics/{0}/{0}.fastq.raw.stats'.format(sample))
        #ls.append('analysis/statistics/{0}/{0}.fastq.clean.stats'.format(sample))

        ls.append('analysis/FastqToSam/{0}/{0}.unmapped.bam'.format(sample))

        ls.append('analysis/MarkIlluminaAdapters/{0}/{0}.MarkIlluminaAdapters.bam'.format(sample))
        ls.append('analysis/MarkIlluminaAdapters/{0}/{0}.MarkIlluminaAdapters.txt'.format(sample))

        ls.append('analysis/SamToFastq/{0}/{0}.interleave.fq'.format(sample))

        ls.append('analysis/bwa_mem/{0}/{0}.mapped.bam'.format(sample))
        ls.append('analysis/bwa_mem/{0}/{0}.mapped.sorted.bam'.format(sample))

        ls.append('analysis/MarkDuplicates/{0}/{0}.dedup.bam'.format(sample))
        ls.append('analysis/MarkDuplicates/{0}/{0}.dedup.sorted.bam'.format(sample))
        ls.append('analysis/MarkDuplicates/{0}/{0}.dedup.txt'.format(sample))
        ls.append('analysis/MarkDuplicates/{0}/{0}.final.bam'.format(sample))

        ls.append('analysis/RecalibrateBaseQualityScores/{0}/{0}.baserecalibrator_data.table'.format(sample))
        ls.append('analysis/RecalibrateBaseQualityScores/{0}/{0}.BQSR.bam'.format(sample))

        ls.append('analysis/HaplotypeCaller/{0}/{0}.g.vcf.gz'.format(sample))
        ls.append('analysis/HaplotypeCaller/{0}/{0}.vcf'.format(sample))
        ls.append('analysis/HaplotypeCaller/{0}/{0}.final.vcf'.format(sample))

        ls.append('analysis/snpeff/{0}/{0}.snpeff.vcf'.format(sample))
        ls.append('analysis/snpeff/{0}/{0}.snpeff.html'.format(sample))

        ls.append("analysis/upload/{0}/{0}.final.bam".format(sample))
        ls.append("analysis/upload/{0}/{0}.g.vcf.gz".format(sample))
        ls.append("analysis/upload/{0}/{0}.final.vcf".format(sample))
        ls.append('analysis/upload/{0}/{0}.snpeff.vcf'.format(sample))
        ls.append('analysis/upload/{0}/{0}.snpeff.html'.format(sample))

    return ls

rule all:
    input: get_outputs(config, samples)

rule prep_ref_genome_bwa_index:
    input:
        genome_fasta = config['reference']['genome_fasta']
    output:
        bwa_index_amb = '{0}.amb'.format(config['reference']['genome_fasta']),
        bwa_index_ann = '{0}.ann'.format(config['reference']['genome_fasta']),
        bwa_index_bwt = '{0}.bwt'.format(config['reference']['genome_fasta']),
        bwa_index_pac = '{0}.pac'.format(config['reference']['genome_fasta']),
        bwa_index_sa = '{0}.sa'.format(config['reference']['genome_fasta'])
    shell:
        "bwa index -a bwtsw {input.genome}"

rule prep_ref_genome_fai:
    input:
        genome_fasta = config['reference']['genome_fasta']
    output:
        genome_fai = '{0}.fai'.format(config['reference']['genome_fasta'])
    shell:
        "samtools faidx {input.genome_fasta}"

rule prep_ref_genome_dict:
    input:
        genome_fasta = config['reference']['genome_fasta']
    output:
        genome_dict = '{0}.dict'.format(config['reference']['genome_fasta'].rstrip('.fasta'))
    shell:
        "picard CreateSequenceDictionary"
        " --REFERENCE {input.genome_fasta}"
        " --OUTPUT {output.genome_dict}"

#rule run_fastqc:
#    input:
#        fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
#        fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
#    output:
#        fastqc_1="analysis/fastqc/{sample}/{sample}_1_fastqc.zip",
#        fastqc_2="analysis/fastqc/{sample}/{sample}_2_fastqc.zip"
#    wildcard_constraints:
#        sample="[^/]+"
#    params:
#        path="analysis/fastqc/{sample}/"
#    threads: 8
#    shell:
#        "fastqc -t {threads}"
#        " --extract"
#        " -o {params.path} {input.fastq_1} {input.fastq_2}"

#rule run_fastq_raw_stats:
#    input:
#        fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
#        fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
#    output:
#        stats="analysis/statistics/{sample}/{sample}.fastq.raw.stats"
#    params:
#        outdir="analysis/statistics/{sample}"
#    wildcard_constraints:
#        sample="[^/]+"
#    shell:
#        "mkdir -p {params.outdir} && "
#        "tool/FastqStatExtractFastqGZ.pairedend.py"
#        " {input.fastq_1} {input.fastq_2} {output.stats}"

rule run_sickle:
    input:
        fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
        fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
    output:
        fastq_1="analysis/sickle/{sample}/{sample}_R1.fastq.gz",
        fastq_2="analysis/sickle/{sample}/{sample}_R2.fastq.gz",
        fastq_s="analysis/sickle/{sample}/{sample}_se.fastq.gz",
    wildcard_constraints:
        sample="[^/]+"
    params:
        log_err="analysis/sickle/{sample}/{sample}.log.err",
        log_out="analysis/sickle/{sample}/{sample}.log.out",
    shell:
        "sickle pe"
        " -f {input.fastq_1} -r {input.fastq_2}"
        " -t sanger -q 20"
        " -l 20 "
        " -g"
        " -o {output.fastq_1} -p {output.fastq_2}"
        " -s {output.fastq_s}"
        " 2> {params.log_err} 1> {params.log_out}"

#rule run_fastq_clean_stats:
#    input:
#        fastq_1="analysis/sickle/{sample}/{sample}_R1.fastq.gz",
#        fastq_2="analysis/sickle/{sample}/{sample}_R2.fastq.gz",
#    output:
#        stats="analysis/statistics/{sample}/{sample}.fastq.clean.stats"
#    params:
#        outdir="analysis/statistics/{sample}"
#    wildcard_constraints:
#        sample="[^/]+"
#    shell:
#        "mkdir -p {params.outdir} && "
#        "tool/FastqStatExtractFastqGZ.pairedend.py"
#        " {input.fastq_1} {input.fastq_2} {output.stats}"

rule run_FastqToSam:
    input:
        fastq_1="analysis/sickle/{sample}/{sample}_R1.fastq.gz",
        fastq_2="analysis/sickle/{sample}/{sample}_R2.fastq.gz",
    output:
        unmapped_bam="analysis/FastqToSam/{sample}/{sample}.unmapped.bam"
    wildcard_constraints:
        sample="[^/]+"
    params:
        sample="{sample}",
        tmp_dir="analysis/FastqToSam/{sample}/tmp",
        log_err="analysis/FastqToSam/{sample}/{sample}.log.err",
        log_out="analysis/FastqToSam/{sample}/{sample}.log.out",
    shell:
        "picard FastqToSam"
        " FASTQ={input.fastq_1}"
        " FASTQ2={input.fastq_2}"
        " OUTPUT={output.unmapped_bam}"
        " READ_GROUP_NAME={params.sample}"
        " SAMPLE_NAME={params.sample}"
        " LIBRARY_NAME={params.sample}"
        #" PLATFORM_UNIT="
        " PLATFORM=ILLUMINA"
        " SEQUENCING_CENTER=TheragenBio"
        " TMP_DIR={params.tmp_dir}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_MarkIlluminaAdapters:
    input:
        unmapped_bam="analysis/FastqToSam/{sample}/{sample}.unmapped.bam"
    output:
        marked_bam="analysis/MarkIlluminaAdapters/{sample}/{sample}.MarkIlluminaAdapters.bam",
        marked_txt="analysis/MarkIlluminaAdapters/{sample}/{sample}.MarkIlluminaAdapters.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        tmp_dir="analysis/MarkIlluminaAdapters/{sample}/tmp",
        log_err="analysis/MarkIlluminaAdapters/{sample}/{sample}.log.err",
        log_out="analysis/MarkIlluminaAdapters/{sample}/{sample}.log.out",
    shell:
        "picard MarkIlluminaAdapters"
        " I={input.unmapped_bam}"
        " O={output.marked_bam}"
        " M={output.marked_txt}"
        " TMP_DIR={params.tmp_dir}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_SamToFastq:
    input:
        marked_bam="analysis/MarkIlluminaAdapters/{sample}/{sample}.MarkIlluminaAdapters.bam",
    output:
        interleave_fq="analysis/SamToFastq/{sample}/{sample}.interleave.fq"
    wildcard_constraints:
        sample="[^/]+"
    params:
        tmp_dir="analysis/SamToFastq/{sample}/tmp",
        log_err="analysis/SamToFastq/{sample}/{sample}.log.err",
        log_out="analysis/SamToFastq/{sample}/{sample}.log.out",
    shell:
        "picard SamToFastq"
        " I={input.marked_bam}"
        " FASTQ={output.interleave_fq}"
        " CLIPPING_ATTRIBUTE=XT"
        " CLIPPING_ACTION=2"
        " INTERLEAVE=true"
        " NON_PF=true"
        " TMP_DIR={params.tmp_dir}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_bwa_mem:
    input:
        interleave_fq="analysis/SamToFastq/{sample}/{sample}.interleave.fq"
    output:
        mapped_bam="analysis/bwa_mem/{sample}/{sample}.mapped.bam",
    wildcard_constraints:
        sample="[^/]+"
    threads: 8
    params:
        idxbase=config["reference"]["genome_fasta"],
        log_err="analysis/bwa_mem/{sample}/{sample}.log.err",
        log_out="analysis/bwa_mem/{sample}/{sample}.log.out",
    shell:
        "bwa mem"
        " -M -t {threads}"
        " -p"
        " -o {output.mapped_bam}"
        " {params.idxbase}"
        " {input.interleave_fq}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_SortSam_mapped:
    input:
        mapped_bam="analysis/bwa_mem/{sample}/{sample}.mapped.bam",
    output:
        sorted_bam="analysis/bwa_mem/{sample}/{sample}.mapped.sorted.bam",
    wildcard_constraints:
        sample="[^/]+"
    params:
        log_err="analysis/bwa_mem/{sample}/{sample}.SortSam.log.err",
        log_out="analysis/bwa_mem/{sample}/{sample}.SortSam.log.out",
    shell:
        "picard SortSam"
        " I={input.mapped_bam}"
        " O={output.sorted_bam}"
        " SORT_ORDER=coordinate"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_MarkDuplicates:
    input:
        mapped_bam="analysis/bwa_mem/{sample}/{sample}.mapped.sorted.bam",
    output:
        dedup_bam="analysis/MarkDuplicates/{sample}/{sample}.dedup.bam",
        dedup_txt="analysis/MarkDuplicates/{sample}/{sample}.dedup.txt",
    wildcard_constraints:
        sample="[^/]+"
    params:
        tmp_dir="analysis/MarkDuplicates/{sample}/tmp",
        log_err="analysis/MarkDuplicates/{sample}/{sample}.log.err",
        log_out="analysis/MarkDuplicates/{sample}/{sample}.log.out",
    shell:
        "picard MarkDuplicates"
        " -Xmx50g"
        " INPUT={input.mapped_bam}"
        " OUTPUT={output.dedup_bam}"
        " METRICS_FILE={output.dedup_txt}"
        " OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
        " CREATE_INDEX=true"
        " CREATE_MD5_FILE=true"
        " REMOVE_SEQUENCING_DUPLICATES=true"
        " TMP_DIR={params.tmp_dir}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_SortSam_dedup:
    input:
        dedup_bam="analysis/MarkDuplicates/{sample}/{sample}.dedup.bam",
    output:
        sorted_bam="analysis/MarkDuplicates/{sample}/{sample}.dedup.sorted.bam",
    wildcard_constraints:
        sample="[^/]+"
    params:
        log_err="analysis/MarkDuplicates/{sample}/{sample}.SortSam.log.err",
        log_out="analysis/MarkDuplicates/{sample}/{sample}.SortSam.log.out",
    shell:
        "picard SortSam"
        " I={input.dedup_bam}"
        " O={output.sorted_bam}"
        " SORT_ORDER=coordinate"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_AddOrReplaceReadGroups:
    input:
        sorted_bam="analysis/MarkDuplicates/{sample}/{sample}.dedup.sorted.bam",
    output:
        final_bam="analysis/MarkDuplicates/{sample}/{sample}.final.bam",
    wildcard_constraints:
        sample="[^/]+"
    params:
        sample_name="{sample}",
        log_err="analysis/MarkDuplicates/{sample}/{sample}.AddOrReplaceReadGroups.log.err",
        log_out="analysis/MarkDuplicates/{sample}/{sample}.AddOrReplaceReadGroups.log.out",
    shell:
        "picard AddOrReplaceReadGroups"
        " --INPUT {input.sorted_bam}"
        " --OUTPUT {output.final_bam}"
        " --RGLB lib1"
        " --RGPL ILLUMINA"
        " --RGPU unit1"
        " --RGSM {params.sample_name}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_BaseRecalibrator:
    input:
        final_bam="analysis/MarkDuplicates/{sample}/{sample}.final.bam",
        genome_dict = '{0}.dict'.format(config['reference']['genome_fasta'].rstrip('.fasta'))
    output:
        baserecal="analysis/RecalibrateBaseQualityScores/{sample}/{sample}.baserecalibrator_data.table"
    wildcard_constraints:
        sample="[^/]+"
    params:
        log_err="analysis/RecalibrateBaseQualityScores/{sample}/{sample}.BaseRecalibrator.log.err",
        log_out="analysis/RecalibrateBaseQualityScores/{sample}/{sample}.BaseRecalibrator.log.out",
    shell:
        "gatk BaseRecalibrator"
        " -I {input.final_bam}"
        " -R {config[reference][genome_fasta]}"
        " --known-sites {config[reference][known_site_1]}"
        " --known-sites {config[reference][known_site_2]}"
        " -O {output.baserecal}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_ApplyBQSR:
    input:
        final_bam="analysis/MarkDuplicates/{sample}/{sample}.final.bam",
        baserecal="analysis/RecalibrateBaseQualityScores/{sample}/{sample}.baserecalibrator_data.table"
    output:
        bqsr_bam="analysis/RecalibrateBaseQualityScores/{sample}/{sample}.BQSR.bam"
    wildcard_constraints:
        sample="[^/]+"
    params:
        log_err="analysis/RecalibrateBaseQualityScores/{sample}/{sample}.ApplyBQSR.log.err",
        log_out="analysis/RecalibrateBaseQualityScores/{sample}/{sample}.ApplyBQSR.log.out",
    shell:
        "gatk ApplyBQSR"
        " -I {input.final_bam}"
        " -R {config[reference][genome_fasta]}"
        " --bqsr-recal-file {input.baserecal}"
        " -O {output.bqsr_bam}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_HaplotypeCaller:
    input:
        bqsr_bam="analysis/RecalibrateBaseQualityScores/{sample}/{sample}.BQSR.bam"
    output:
        gvcf="analysis/HaplotypeCaller/{sample}/{sample}.g.vcf.gz"
    wildcard_constraints:
        sample="[^/]+"
    threads: 8
    params:
        log_err="analysis/HaplotypeCaller/{sample}/{sample}.log.err",
        log_out="analysis/HaplotypeCaller/{sample}/{sample}.log.out",
    shell:
        "gatk HaplotypeCaller"
        " -I {input.bqsr_bam}"
        " -R {config[reference][genome_fasta]}"
        " -O {output.gvcf}"
        " -ERC GVCF"
        " --native-pair-hmm-threads {threads}"
        " 2> {params.log_err} 1> {params.log_out}"

rule run_HaplotypeCaller_VCF:
    input:
        gvcf="analysis/HaplotypeCaller/{sample}/{sample}.g.vcf.gz"
        #gvcf="analysis/gatk4/{sample}/{sample}.recal.haplotypecall.g.vcf",
    output:
        vcf="analysis/HaplotypeCaller/{sample}/{sample}.vcf",
        tmp=temp("analysis/HaplotypeCaller/{sample}/temp"),
        #vcf="analysis/gatk4/{sample}/{sample}.recal.haplotypecall.vcf",
        #tmp=temp("analysis/gatk4/{sample}/temp"),
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "mkdir -p {output.tmp} && "
        "gatk GenotypeGVCFs"
        " --tmp-dir {output.tmp}"
        " -R {config[reference][genome_fasta]}"
        #" --intervals {config[variant_target_bed]}"
        " -V {input.gvcf} "
        " -O {output.vcf}"

rule gatk_filter:
    input:
        vcf="analysis/HaplotypeCaller/{sample}/{sample}.vcf",
        #vcf="analysis/gatk4/{sample}/{sample}.recal.haplotypecall.vcf",
    output:
        tmp=temp("analysis/HaplotypeCaller/{sample}/temp"),
        filtvcf="analysis/HaplotypeCaller/{sample}/{sample}.final.vcf",
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "mkdir -p {output.tmp} && "
        "gatk VariantFiltration"
        " --tmp-dir {output.tmp}"
        " -V {input.vcf}"
        " -O {output.filtvcf}"

rule AnnotateVCF:
    input:
        vcf="analysis/HaplotypeCaller/{sample}/{sample}.final.vcf"
    output:
        vcf="analysis/snpeff/{sample}/{sample}.snpeff.vcf",
        html="analysis/snpeff/{sample}/{sample}.snpeff.html"
    wildcard_constraints:
        sample="[^/]+"
    params:
        updownstream=2000,
        sampleid="{sample}"
    shell:
        "snpEff -geneId"
        " -noInteraction"
        " -noMotif"
        " -noNextProt"
        " -no-intergenic"
        " -ud {params.updownstream}"
        " -canon"
        " -c {config[snpeff_config]}"
        " -v {config[snpeff_db]}"
        " -s {output.html}"
        " -o vcf"
        " {input.vcf} | "
        "SnpSift varType - | "
        "SnpSift annotate -a -ID  {config[dbsnp_db]} - | "
        "SnpSift annotate -a -noID  {config[clinvar_db]} - "
        " > {output.vcf}"
        #"SnpSift DbNsfp -db {config[dbnsftp_db]} -f {config[snpeff_annot_columns]} -  | "
        #"SnpSift annotate -a -noID  {config[cosmic_db]} - | "
        #"SnpSift annotate -a -noID  {config[gnomad_db]} - | "
        #"SnpSift annotate -a -noID  {config[krgdb]} - "

rule symlink_upload:
    input:
        final_bam="analysis/MarkDuplicates/{sample}/{sample}.final.bam",
        gvcf="analysis/HaplotypeCaller/{sample}/{sample}.g.vcf.gz",
        vcf="analysis/HaplotypeCaller/{sample}/{sample}.final.vcf",
        snpeff_vcf="analysis/snpeff/{sample}/{sample}.snpeff.vcf",
        snpeff_html="analysis/snpeff/{sample}/{sample}.snpeff.html",
    output:
        final_bam="analysis/upload/{sample}/{sample}.final.bam",
        gvcf="analysis/upload/{sample}/{sample}.g.vcf.gz",
        vcf="analysis/upload/{sample}/{sample}.final.vcf",
        snpeff_vcf="analysis/upload/{sample}/{sample}.snpeff.vcf",
        snpeff_html="analysis/upload/{sample}/{sample}.snpeff.html",
    params:
        wkdir=os.path.abspath(config["workdir"])
    shell:
        "ln -s {params.wkdir}/{input.final_bam} {params.wkdir}/{output.final_bam} && "
        "ln -s {params.wkdir}/{input.gvcf} {params.wkdir}/{output.gvcf} && "
        "ln -s {params.wkdir}/{input.vcf} {params.wkdir}/{output.vcf} && "
        "ln -s {params.wkdir}/{input.snpeff_vcf} {params.wkdir}/{output.snpeff_vcf} && "
        "ln -s {params.wkdir}/{input.snpeff_html} {params.wkdir}/{output.snpeff_html}"


