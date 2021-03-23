

#import yaml


class Config:
    def __init__(self, rawdata_path, samples, outfn):
        self.rawdata_path = rawdata_path
        self.samples = samples
        self.outfn = outfn

    def make_yaml(self):
        outfh = open(self.outfn, 'w')
        print('', file=outfh)
        #print('workdir: /BiO/BioPeople/siyoo/de_novo_mutation_calling_pipeline', file=outfh)
        print('workdir: /data02/project/TBD200742_WGRS_Human_siyoo_202011/de_novo_mutation_calling_pipeline', file=outfh)
        print('reference:', file=outfh)
        print('  genome_fasta: reference/Homo_sapiens_assembly38.fasta', file=outfh)
        print('  known_site_1: reference/Homo_sapiens_assembly38.dbsnp138.vcf', file=outfh)
        print('  known_site_2: reference/Mills_and_1000G_gold_standard.indels.hg38.vcf', file=outfh)
        print('snpeff_config: /TBI/People/tbi/siyoo/miniconda3/envs/de_novo_mutation_calling_pipeline/share/snpeff-5.0-0/snpEff.config', file=outfh)
        print('snpeff_db: hg38', file=outfh)
        #print('dbnsftp_db: /TBI/Share/BioPeople/shsong/BioResources/DBs/dbNSFP/dbNSFP4.0b2a.hg19.txt.gz', file=outfh)
        #print('cosmic_db: /TBI/Share/BioPeople/shsong/BioResources/DBs/cosmic/v79/VCF/CosmicCodingMuts.anno.vcf.gz', file=outfh)
        print('dbsnp_db: /TBI/Share/BioPeople/shsong/BioResources/DBs/dbsnp/v150/GRCh38/00-All.vcf.gz', file=outfh)
        print('clinvar_db: /TBI/Share/BioPeople/shsong/BioResources/DBs/clinvar/20200905/GRCh38/clinvar.vcf.gz', file=outfh)
        #print('gnomad_db: /TBI/Share/BioPeople/shsong/BioResources/DBs/gnomad/r2.1.1/gnomad.exomes.r2.1.1.sites.filter.20190612.vcf.gz', file=outfh)
        #print('krgdb: /TBI/Share/BioPeople/shsong/BioResources/DBs/KRGDB/KRGDB_COMMON_SNP_INDEL.vcf.gz', file=outfh)
        #print('snpeff_annot_columns: "chr,ref,alt,aaref,aaalt,rs_dbSNP151,aapos,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,Uniprot_acc,Uniprot_entry,APPRIS,GENCODE_basic,TSL,VEP_canonical,cds_strand,refcodon,codonpos,codon_degeneracy,Ancestral_allele,AltaiNeandertal,Denisova,VindijiaNeandertal,SIFT_score,SIFT_converted_rankscore,SIFT_pred,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,MVP_score,MVP_rankscore,MPC_score,MPC_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,DEOGEN2_score,DEOGEN2_rankscore,DEOGEN2_pred,Aloft_Fraction_transcripts_affected,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,CADD_raw,CADD_raw_rankscore,CADD_phred,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,fathmm-XF_coding_score,fathmm-XF_coding_rankscore,fathmm-XF_coding_pred,Eigen-raw_coding,Eigen-raw_coding_rankscore,Eigen-pred_coding,Eigen-PC-raw_coding,Eigen-PC-raw_coding_rankscore,Eigen-PC-phred_coding,GenoCanyon_score,GenoCanyon_rankscore,integrated_fitCons_score,integrated_fitCons_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_rankscore,HUVEC_confidence_value,LINSIGHT,LINSIGHT_rankscore,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian,phyloP30way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons30way_mammalian,phastCons30way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,bStatistic,bStatistic_rankscore,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,UK10K_AC,UK10K_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,ExAC_nonTCGA_AC,ExAC_nonTCGA_AF,ExAC_nonTCGA_Adj_AC,ExAC_nonTCGA_Adj_AF,ExAC_nonTCGA_AFR_AC,ExAC_nonTCGA_AFR_AF,ExAC_nonTCGA_AMR_AC,ExAC_nonTCGA_AMR_AF,ExAC_nonTCGA_EAS_AC,ExAC_nonTCGA_EAS_AF,ExAC_nonTCGA_FIN_AC,ExAC_nonTCGA_FIN_AF,ExAC_nonTCGA_NFE_AC,ExAC_nonTCGA_NFE_AF,ExAC_nonTCGA_SAS_AC,ExAC_nonTCGA_SAS_AF,ExAC_nonpsych_AC,ExAC_nonpsych_AF,ExAC_nonpsych_Adj_AC,ExAC_nonpsych_Adj_AF,ExAC_nonpsych_AFR_AC,ExAC_nonpsych_AFR_AF,ExAC_nonpsych_AMR_AC,ExAC_nonpsych_AMR_AF,ExAC_nonpsych_EAS_AC,ExAC_nonpsych_EAS_AF,ExAC_nonpsych_FIN_AC,ExAC_nonpsych_FIN_AF,ExAC_nonpsych_NFE_AC,ExAC_nonpsych_NFE_AF,ExAC_nonpsych_SAS_AC,ExAC_nonpsych_SAS_AF,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,Interpro_domain,GTEx_V7_gene,GTEx_V7_tissue,Geuvadis_eQTL_target_gene"', file=outfh)
        print('', file=outfh)
        print('samples:', file=outfh)
        for sample in self.samples:
            print('  {0}:'.format(sample), file=outfh)
            print('    fastq_1: {0}/{1}_1.fq.gz'.format(self.rawdata_path, sample), file=outfh)
            print('    fastq_2: {0}/{1}_2.fq.gz'.format(self.rawdata_path, sample), file=outfh)
        outfh.close()

def main():
    #samples = ['sample_K-ABC-00002', 'sample_K-ABC-00011', 'sample_K-ABC-00015', 'sample_K-ABC-00016',
    #           'sample_K-ABC-00017', 'sample_K-ABC-00020', 'sample_K-ABC-00022', 'sample_K-ABC-00026',
    #           'sample_K-ABC-00027', 'sample_K-ABC-00034', 'sample_K-ABC-00035', 'sample_K-ABC-00036',
    #           'sample_K-ABC-00037', 'sample_K-ABC-00038', 'sample_K-ABC-00040', 'sample_K-ABC-00044',
    #           'sample_K-ABC-00087']
    #samples = ['sample_K-CTC-00001', 'sample_K-CTC-00002', 'sample_K-CTC-00003', 'sample_K-CTC-00004',
    #           'sample_K-CTC-00005', 'sample_K-CTC-00006', 'sample_K-CTC-00007', 'sample_K-CTC-00008',
    #           'sample_K-CTC-00009', 'sample_K-CTC-00010', 'sample_K-CTC-00011', 'sample_K-CTC-00012',
    #           'sample_K-CTC-00013', 'sample_K-CTC-00014', 'sample_K-CTC-00015', 'sample_K-CTC-00016',
    #           'sample_K-CTC-00017', 'sample_K-CTC-00018', 'sample_K-CTC-00019', 'sample_K-CTC-00020']
    #rawdata_path = 'rawdata/TBD200889_10869_2531_20201210'
    #samples = ['sample_K-CTC-00021', 'sample_K-CTC-00022', 'sample_K-CTC-00023', 'sample_K-CTC-00024',
    #           'sample_K-CTC-00025', 'sample_K-CTC-00026', 'sample_K-CTC-00027', 'sample_K-CTC-00028',
    #           'sample_K-CTC-00029', 'sample_K-CTC-00030', 'sample_K-CTC-00031', 'sample_K-CTC-00032',
    #           'sample_K-CTC-00033', 'sample_K-CTC-00034', 'sample_K-CTC-00035', 'sample_K-CTC-00036',
    #           'sample_K-CTC-00037', 'sample_K-CTC-00038', 'sample_K-CTC-00039']
    #rawdata_path = 'rawdata/TBD201002_11148_2778_20210112'
    #
    rawdata_path = 'rawdata/TBD210101_11481_20210303'
    samples = ['sample_K-ABC-00209', 'sample_K-ABC-00278']

    outfn = 'config.yaml'
    config = Config(rawdata_path, samples, outfn)
    config.make_yaml()

if __name__=='__main__':
    main()
