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
    ls.append('{0}.dict'.format(config['reference']['genome_fasta']))

    for sample in samples:
        ls.append('analysis/fastqc/{0}/{0}_1_fastqc.zip'.format(sample))
        ls.append('analysis/fastqc/{0}/{0}_2_fastqc.zip'.format(sample))

        ls.append('analysis/sickle/{0}/{0}_R1.fastq.gz'.format(sample))
        ls.append('analysis/sickle/{0}/{0}_R2.fastq.gz'.format(sample))
        ls.append('analysis/sickle/{0}/{0}_se.fastq.gz'.format(sample))

        ls.append('analysis/statistics/{0}/{0}.fastq.raw.stats'.format(sample))
        ls.append('analysis/statistics/{0}/{0}.fastq.clean.stats'.format(sample))

        ls.append('analysis/FastqToSam/{0}/{0}.unmapped.bam'.format(sample))

        ls.append('analysis/MarkIlluminaAdapters/{0}/{0}.MarkIlluminaAdapters.bam'.format(sample))
        ls.append('analysis/MarkIlluminaAdapters/{0}/{0}.MarkIlluminaAdapters.txt'.format(sample))

        ls.append('analysis/SamToFastq/{0}/{0}.interleave.fq'.format(sample))

        ls.append('analysis/bwa_mem/{0}/{0}.mapped.bam'.format(sample))

        ls.append('analysis/MarkDuplicates/{0}/{0}.dedup.bam'.format(sample))
        ls.append('analysis/MarkDuplicates/{0}/{0}.dedup.txt'.format(sample))

        ls.append('analysis/SortSam/{0}/{0}.sorted.bam'.format(sample))

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
        genome_dict = '{0}.dict'.format(config['reference']['genome_fasta'])
    shell:
        "picard CreateSequenceDictionary"
        " --REFERENCE {input.genome_fasta}"
        " --OUTPUT {output.genome_dict}"

rule run_fastqc:
    input:
        fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
        fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
    output:
        fastqc_1="analysis/fastqc/{sample}/{sample}_1_fastqc.zip",
        fastqc_2="analysis/fastqc/{sample}/{sample}_2_fastqc.zip"
    wildcard_constraints:
        sample="[^/]+"
    params:
        path="analysis/fastqc/{sample}/"
    threads: 8
    shell:
        "fastqc -t {threads}"
        " --extract"
        " -o {params.path} {input.fastq_1} {input.fastq_2}"

rule run_fastq_raw_stats:
    input:
        fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
        fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
    output:
        stats="analysis/statistics/{sample}/{sample}.fastq.raw.stats"
    params:
        outdir="analysis/statistics/{sample}"
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "mkdir -p {params.outdir} && "
        "tool/FastqStatExtractFastqGZ.pairedend.py"
        " {input.fastq_1} {input.fastq_2} {output.stats}"

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
    shell:
        "sickle pe"
        " -f {input.fastq_1} -r {input.fastq_2}"
        " -t sanger -q 20"
        " -l 20 "
        " -g"
        " -o {output.fastq_1} -p {output.fastq_2}"
        " -s {output.fastq_s}"

rule run_fastq_clean_stats:
    input:
        fastq_1="analysis/sickle/{sample}/{sample}_R1.fastq.gz",
        fastq_2="analysis/sickle/{sample}/{sample}_R2.fastq.gz",
    output:
        stats="analysis/statistics/{sample}/{sample}.fastq.clean.stats"
    params:
        outdir="analysis/statistics/{sample}"
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "mkdir -p {params.outdir} && "
        "tool/FastqStatExtractFastqGZ.pairedend.py"
        " {input.fastq_1} {input.fastq_2} {output.stats}"

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

rule run_MarkIlluminaAdapters:
    input:
        unmapped_bam="analysis/FastqToSam/{sample}/{sample}.unmapped.bam"
    output:
        marked_bam="analysis/MarkIlluminaAdapters/{sample}/{sample}.MarkIlluminaAdapters.bam",
        marked_txt="analysis/MarkIlluminaAdapters/{sample}/{sample}.MarkIlluminaAdapters.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        tmp_dir="analysis/MarkIlluminaAdapters/{sample}/tmp"
    shell:
        "picard MarkIlluminaAdapters"
        " I={input.unmapped_bam}"
        " O={output.marked_bam}"
        " M={output.marked_txt}"
        " TMP_DIR={params.tmp_dir}"

rule run_SamToFastq:
    input:
        marked_bam="analysis/MarkIlluminaAdapters/{sample}/{sample}.MarkIlluminaAdapters.bam",
    output:
        interleave_fq="analysis/SamToFastq/{sample}/{sample}.interleave.fq"
    wildcard_constraints:
        sample="[^/]+"
    params:
        tmp_dir="analysis/SamToFastq/{sample}/tmp"
    shell:
        "picard SamToFastq"
        " I={input.marked_bam}"
        " FASTQ={output.interleave_fq}"
        " CLIPPING_ATTRIBUTE=XT"
        " CLIPPING_ACTION=2"
        " INTERLEAVE=true"
        " NON_PF=true"
        " TMP_DIR={params.tmp_dir}"

rule run_bwa_mem:
    input:
        interleave_fq="analysis/SamToFastq/{sample}/{sample}.interleave.fq"
    output:
        mapped_bam="analysis/bwa_mem/{sample}/{sample}.mapped.bam",
    wildcard_constraints:
        sample="[^/]+"
    threads: 8
    params:
        idxbase=config["reference"]["genome_fasta"]
    shell:
        "bwa mem"
        " -M -t {threads}"
        " -p"
        " -o {output.mapped_bam}"
        " {params.idxbase}"
        " {input.interleave_fq}"

rule run_MarkDuplicates:
    input:
        mapped_bam="analysis/bwa_mem/{sample}/{sample}.mapped.bam",
    output:
        dedup_bam="analysis/MarkDuplicates/{sample}/{sample}.dedup.bam",
        dedup_txt="analysis/MarkDuplicates/{sample}/{sample}.dedup.txt",
    wildcard_constraints:
        sample="[^/]+"
    params:
        tmp_dir="analysis/MarkDuplicates/{sample}/tmp"
    shell:
        "picard MarkDuplicates"
        " INPUT={input.mapped_bam}"
        " OUTPUT={output.dedup_bam}"
        " METRICS_FILE={output.dedup_txt}"
        " OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
        " CREATE_INDEX=true"
        " ASSUME_SORT_ORDER=unsorted"
        " TMP_DIR={params.tmp_dir}"

rule run_SortSam:
    input:
        dedup_bam="analysis/MarkDuplicates/{sample}/{sample}.dedup.bam",
    output:
        sorted_bam="analysis/SortSam/{sample}/{sample}.sorted.bam",
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "picard SortSam"
        " I={input.dedup_bam}"
        " O={output.sorted_bam}"
        " SORT_ORDER=coordinate"



