#!/usr/bin/python

import os
import yaml

configfile: "config.yaml"

wkdir = "/BiO/BioPeople/siyoo/de_novo_mutation_calling_pipeline"
samples = ['sample1','sample2','sample3']

def get_outputs(config, wkdir, samples):
    ls = list()
    ls.append('{0}.amb'.format(config['reference']['genome_fasta']))
    ls.append('{0}.ann'.format(config['reference']['genome_fasta']))
    ls.append('{0}.bwt'.format(config['reference']['genome_fasta']))
    ls.append('{0}.pac'.format(config['reference']['genome_fasta']))
    ls.append('{0}.sa'.format(config['reference']['genome_fasta']))
    ls.append('{0}.fai'.format(config['reference']['genome_fasta']))
    ls.append('{0}.dict'.format(config['reference']['genome_fasta']))

    for sample in samples:
        ls.append('analysis/fastqc/{0}/{0}_R1_fastqc.zip'.format(sample))
        ls.append('analysis/fastqc/{0}/{0}_R2_fastqc.zip'.format(sample))

        ls.append('analysis/sickle/{0}/{0}_R1.fastq.gz'.format(sample))
        ls.append('analysis/sickle/{0}/{0}_R2.fastq.gz'.format(sample))
        ls.append('analysis/sickle/{0}/{0}_se.fastq.gz'.format(sample))

        ls.append('analysis/statistics/{0}/{0}.fastq.raw.stats'.format(sample))
        ls.append('analysis/statistics/{0}/{0}.fastq.clean.stats'.format(sample))

    return ls

rule all:
    input: get_outputs(config, wkdir, samples)

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
        fastqc_1="analysis/fastqc/{sample}/{sample}_R1_fastqc.zip",
        fastqc_2="analysis/fastqc/{sample}/{sample}_R2_fastqc.zip"
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
        " -I 20 "
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





