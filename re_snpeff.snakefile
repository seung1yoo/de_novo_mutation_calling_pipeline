

workdir: '/data02/project/TBD200742_WGRS_Human_siyoo_202011/de_novo_mutation_calling_pipeline/analysis'

samples=["sample_K-CTC-00021", "sample_K-CTC-00022", "sample_K-CTC-00023", "sample_K-CTC-00024", "sample_K-CTC-00025",
         "sample_K-CTC-00026", "sample_K-CTC-00027", "sample_K-CTC-00028", "sample_K-CTC-00029", "sample_K-CTC-00030",
         "sample_K-CTC-00031", "sample_K-CTC-00032", "sample_K-CTC-00033", "sample_K-CTC-00034", "sample_K-CTC-00035",
         "sample_K-CTC-00036", "sample_K-CTC-00037", "sample_K-CTC-00038", "sample_K-CTC-00039"]

def get_outputs(samples):
    ls = list()
    for sample in samples:
        ls.append(f"snpeff/{sample}/{sample}.snpeff.vcf")
        ls.append(f"snpeff/{sample}/{sample}.snpeff.html")
    return ls


rule all:
    input: get_outputs(samples)

rule AnnotateVCF:
    input:
        vcf="HaplotypeCaller/{sample}/{sample}.final.vcf"
    output:
        vcf="snpeff/{sample}/{sample}.snpeff.vcf",
        html="snpeff/{sample}/{sample}.snpeff.html"
    wildcard_constraints:
        sample="[^/]+"
    params:
        updownstream=2000,
        sampleid="{sample}",
        snpeff_config='/TBI/People/tbi/siyoo/miniconda3/envs/de_novo_mutation_calling_pipeline/share/snpeff-5.0-0/snpEff.config',
        snpeff_db="hg38",
        dbsnp_db="/TBI/Share/BioPeople/shsong/BioResources/DBs/dbsnp/v150/GRCh38/00-All.vcf.gz",
        clinvar_db="/TBI/Share/BioPeople/shsong/BioResources/DBs/clinvar/20200905/GRCh38/clinvar.vcf.gz"
    shell:
        "snpEff -geneId"
        " -noInteraction"
        " -noMotif"
        " -noNextProt"
        " -no-intergenic"
        " -ud {params.updownstream}"
        " -canon"
        " -c {params.snpeff_config}"
        " -v {params.snpeff_db}"
        " -s {output.html}"
        " -o vcf"
        " {input.vcf} | "
        "SnpSift varType - | "
        "SnpSift annotate -a -ID  {params.dbsnp_db} - | "
        "SnpSift annotate -a -noID  {params.clinvar_db} - "
        " > {output.vcf}"
