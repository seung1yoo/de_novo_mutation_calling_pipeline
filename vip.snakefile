#!/usr/bin/python

import os
import yaml

configfile: "config.yaml"

wkdir = "/BiO/BioPeople/siyoo/de_novo_mutation_calling_pipeline"
samples = ['AKTPF_hez_Acv1','AKTPF_hez_Acv2','AKTPF_hez_Acv3',
           'AKTPF_hez_control1','AKTPF_hez_control2','AKTPF_hez_control3',
           'AKTPF_LOH_Acv1','AKTPF_LOH_Acv2','AKTPF_LOH_Acv3',
           'AKTPF_LOH_control1','AKTPF_LOH_control2','AKTPF_LOH_control3']
comps = ['DEG01','DEG02','DEG03','DEG04']

def get_outputs(wkdir, samples, comps):
    ls = list()
    ls.append(os.path.join(wkdir, 'star-genome', 'SAindex'))
    for sample in samples:
        ls.append('analysis/cutadapt/{0}.clean.fq.gz'.format(sample))
        ls.append('analysis/fastqc/{0}_fastqc.zip'.format(sample))
        ls.append('analysis/star/{0}/Aligned.sortedByCoord.out.bam'.format(sample))
        ls.append('analysis/star/{0}/Aligned.sortedByCoord.out.bam.stats'.format(sample))
        ls.append('analysis/star/{0}/Aligned.sortedByCoord.out.bam.bai'.format(sample))
        ls.append('analysis/rseqc/{0}/{0}.DupRate_plot.pdf'.format(sample))
        ls.append('analysis/rseqc/{0}/{0}.DupRate_plot.r'.format(sample))
        ls.append('analysis/rseqc/{0}/{0}.pos.DupRate.xls'.format(sample))
        ls.append('analysis/rseqc/{0}/{0}.seq.DupRate.xls'.format(sample))
        ls.append('analysis/cuffquant/{0}/abundances.cxb'.format(sample))
        ls.append('analysis/cufflinks/{0}/genes.fpkm_tracking'.format(sample))
    for comp in comps:
        ls.append('analysis/cuffdiff/{0}/gene_exp.diff'.format(comp))
    ls.append("analysis/cuffnorm/genes.fpkm_table")
    return ls

def get_cxb(config, comparison, group):
    cxb_s = list()
    for sample in config['comps'][comparison][group].split(','):
        cxb_s.append('analysis/cuffquant/{0}/abundances.cxb'.format(sample))
    return ','.join(cxb_s)

rule all:
    input: get_outputs(wkdir, samples, comps)

rule starGenomeGenerate:
    input:
        genome  = "genome.fa",
        geneset = "genes.gtf"
    output:
        starindex = os.path.join(wkdir,'star-genome','SAindex')
    params:
        stargenomedir = os.path.join(wkdir,'star-genome')
    threads: 10
    message: "Run Reference Genome Building with STAR"
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate"
        " --genomeDir {params.stargenomedir}"
        " --genomeFastaFiles {input.genome}"
        " --sjdbGTFfile {input.geneset}"
        " --limitGenomeGenerateRAM 31000000000"

rule run_cutadapt:
    input:
        raw_fq=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"]
    output:
        clean_fq="analysis/cutadapt/{sample}.clean.fq.gz",
        log_fn="analysis/cutadapt/{sample}.cutadapt.log"
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "unset PYTHONPATH && "
        "cutadapt -q 30"
        " -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        " --minimum-length 30"
        " -o {output.clean_fq}"
        " {input.raw_fq} > {output.log_fn}"

rule run_fastqc:
    input:
        raw_fq=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"]
    output:
        fastqc_fn="analysis/fastqc/{sample}_fastqc.zip"
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "fastqc "
        " --extract"
        " -o analysis/fastqc {input.raw_fq}"

rule run_star:
    input:
        fastq_1="analysis/cutadapt/{sample}.clean.fq.gz",
    output:
        bam="analysis/star/{sample}/Aligned.out.bam",
        sortBam="analysis/star/{sample}/Aligned.sortedByCoord.out.bam",
        tranBam="analysis/star/{sample}/Aligned.toTranscriptome.out.bam",
        readsPerGene="analysis/star/{sample}/ReadsPerGene.out.tab",
        logFinal="analysis/star/{sample}/Log.final.out",
        SJout="analysis/star/{sample}/SJ.out.tab",
        unmapRead1="analysis/star/{sample}/Unmapped.out.mate1",
    wildcard_constraints:
        sample="[^/]+"
    params:
        outdir = 'analysis/star/{sample}',
        prefix = 'analysis/star/{sample}/',
        stargenomedir = os.path.join(wkdir,'star-genome'),
        geneset = 'genes.gtf'
    threads: 8
    shell:
        "mkdir -p {params.outdir} && "
        "STAR --runMode alignReads"
        " --runThreadN {threads}"
        " --genomeDir {params.stargenomedir}"
        " --readFilesIn {input.fastq_1}"
        " --readFilesCommand zcat"
        " --outFileNamePrefix {params.prefix}"
        " --outSAMtype BAM Unsorted SortedByCoordinate"
        " --outSAMstrandField intronMotif"
        " --quantMode TranscriptomeSAM GeneCounts"
        " --sjdbGTFfile {params.geneset}"
        " --outReadsUnmapped Fastx"

rule run_star_idx:
    input:
        sortBam="analysis/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        sortBamIdx="analysis/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "samtools index {input.sortBam}"

rule run_samtools_stats:
    input:
        sortBam="analysis/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        stats="analysis/star/{sample}/Aligned.sortedByCoord.out.bam.stats",
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "samtools stats {input.sortBam} > {output.stats}"

rule run_rseqc_readdup:
    input:
        sortBam="analysis/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        duprate_pdf="analysis/rseqc/{sample}/{sample}.DupRate_plot.pdf",
        duprate_r="analysis/rseqc/{sample}/{sample}.DupRate_plot.r",
        duprate_pos="analysis/rseqc/{sample}/{sample}.pos.DupRate.xls",
        duprate_seq="analysis/rseqc/{sample}/{sample}.seq.DupRate.xls",
    wildcard_constraints:
        sample="[^/]+"
    params:
        outprefix="analysis/rseqc/{sample}/{sample}"
    threads: 1
    shell:
        "unset PYTHONPATH && "
        "read_duplication.py "
        " -i {input.sortBam}"
        " -o {params.outprefix}"

rule run_cuffquant:
    input:
        sortBam="analysis/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        cxb="analysis/cuffquant/{sample}/abundances.cxb"
    wildcard_constraints:
        sample="[^/]+"
    params:
        outdir="analysis/cuffquant/{sample}",
        geneset="genes.gtf"
    threads: 8
    shell:
        "mkdir -p {params.outdir} && "
        " cuffquant "
        " --num-threads {threads}"
        " --output-dir {params.outdir}"
        " {params.geneset}"
        " {input.sortBam}"

rule run_cufflinks:
    input:
        cxb="analysis/cuffquant/{sample}/abundances.cxb"
    output:
        fpkm_fn="analysis/cufflinks/{sample}/genes.fpkm_tracking"
    wildcard_constraints:
        sample="[^/]+"
    params:
        outdir="analysis/cufflinks/{sample}",
        geneset="genes.gtf"
    threads: 8
    shell:
        "cufflinks -p {threads}"
        " -G {params.geneset}"
        " -o {params.outdir}"
        " {input.cxb}"

rule run_cuffnorm:
    input:
        cxb_s=expand("analysis/cuffquant/{sample}/abundances.cxb", sample=samples)
    output:
        fpkm_fn="analysis/cuffnorm/genes.fpkm_table"
    wildcard_constraints:
        sample="[^/]+"
    params:
        labels=','.join(samples),
        outdir="analysis/cuffnorm",
        geneset="genes.gtf"
    threads: 8
    shell:
        " cuffnorm "
        " --num-threads {threads}"
        " --output-dir {params.outdir}"
        " --labels {params.labels}"
        " --library-norm-method classic-fpkm"
        " {params.geneset}"
        " {input.cxb_s}"

rule run_cuffdiff:
    input:
        cxb_s=expand("analysis/cuffquant/{sample}/abundances.cxb", sample=samples)
    output:
        diff_fn="analysis/cuffdiff/{comparison}/gene_exp.diff"
    wildcard_constraints:
        sample="[^/]+",
        comparison="[^/]+"
    params:
        ctrl_id=lambda wildcards: config['comps'][wildcards.comparison]['ctrl'],
        case_id=lambda wildcards: config['comps'][wildcards.comparison]['case'],
        ctrl_cxb=lambda wildcards: get_cxb(config, wildcards.comparison, 'ctrl'),
        case_cxb=lambda wildcards: get_cxb(config, wildcards.comparison, 'case'),
        ctrl_label=lambda wildcards: config['comps'][wildcards.comparison]['ctrl_label'],
        case_label=lambda wildcards: config['comps'][wildcards.comparison]['case_label'],
        outdir="analysis/cuffdiff/{comparison}",
        geneset="genes.gtf"
    threads: 8
    shell:
        "cuffdiff --num-threads {threads}"
        " --output-dir {params.outdir}"
        " --labels {params.ctrl_label},{params.case_label}"
        " --library-norm-method classic-fpkm"
        " {params.geneset}"
        " {params.ctrl_cxb}"
        " {params.case_cxb}"











