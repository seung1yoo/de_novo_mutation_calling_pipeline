

#import yaml


class Config:
    def __init__(self, samples, outfn):
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
        print('  known_site_2: Mills_and_1000G_gold_standard.indels.hg38.vcf', file=outfh)
        print('', file=outfh)
        print('samples:', file=outfh)
        for sample in self.samples:
            print('  {0}:'.format(sample), file=outfh)
            print('    fastq_1: rawdata/{0}_1.fq.gz'.format(sample), file=outfh)
            print('    fastq_2: rawdata/{0}_2.fq.gz'.format(sample), file=outfh)
        outfh.close()

def main():
    samples = ['sample_K-ABC-00002',
               'sample_K-ABC-00011',
               'sample_K-ABC-00015',
               'sample_K-ABC-00016',
               'sample_K-ABC-00017',
               'sample_K-ABC-00020',
               'sample_K-ABC-00022',
               'sample_K-ABC-00026',
               'sample_K-ABC-00027',
               'sample_K-ABC-00034',
               'sample_K-ABC-00035',
               'sample_K-ABC-00036',
               'sample_K-ABC-00037',
               'sample_K-ABC-00038',
               'sample_K-ABC-00040',
               'sample_K-ABC-00044',
               'sample_K-ABC-00087']
    outfn = 'config.yaml'
    config = Config(samples, outfn)
    config.make_yaml()

if __name__=='__main__':
    main()
