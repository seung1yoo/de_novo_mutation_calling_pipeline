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

    for sample in samples:
        pass
        #ls.append('analysis/cutadapt/{0}.clean.fq.gz'.format(sample))
        #ls.append('analysis/fastqc/{0}_fastqc.zip'.format(sample))

    return ls

rule all:
    input: get_outputs(config, wkdir, samples)

rule prep_bwa_index:
    input:
        genome = config['reference']['genome_fasta']
    output:
        amb = '{0}.amb'.format(config['reference']['genome_fasta']),
        ann = '{0}.ann'.format(config['reference']['genome_fasta']),
        bwt = '{0}.bwt'.format(config['reference']['genome_fasta']),
        pac = '{0}.pac'.format(config['reference']['genome_fasta']),
        sa = '{0}.sa'.format(config['reference']['genome_fasta'])
    shell:
        "bwa index -a bwtsw {input.genome}"








