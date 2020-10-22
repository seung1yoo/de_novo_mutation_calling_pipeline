

def star_logFinal_parser(fn):
    in_read = 0
    mapped_reads = 0
    mapped_reads_uniquely = 0
    mapped_reads_multiple = 0
    unmapped_reads = 0
    for line in open(fn):
        if not line.strip():
            continue
        items = line.strip().split('\t')
        if items[0].startswith('Number of input reads'):
            in_read = int(items[1])
        elif items[0].startswith('Uniquely mapped reads number'):
            mapped_reads += int(items[1])
            mapped_reads_uniquely = int(items[1])
        elif items[0].startswith('Number of reads mapped to multiple loci'):
            mapped_reads += int(items[1])
            mapped_reads_multiple += int(items[1])
        elif items[0].startswith('Number of reads mapped to too many loci'):
            mapped_reads += int(items[1])
            mapped_reads_multiple += int(items[1])
        elif items[0].startswith('Number of reads unmapped: too many mismatches'):
            unmapped_reads += int(items[1])
        elif items[0].startswith('Number of reads unmapped: too short'):
            unmapped_reads += int(items[1])
        elif items[0].startswith('Number of reads unmapped: other'):
            unmapped_reads += int(items[1])

    if not mapped_reads in [mapped_reads_uniquely+mapped_reads_multiple]:
        print('ERROR: check the star_logFinal_parser')
        print('if not mapped_reads in [mapped_reads_uniquely+mapped_reads_multiple]:')
        print('in_read', in_read)
        print('mapped_reads', mapped_reads)
        print('mapped_reads_uniquely', mapped_reads_uniquely)
        print('mapped_reads_multiple', mapped_reads_multiple)
        print('unmapped_reads', unmapped_reads)
        sys.exit()
    if not in_read in [mapped_reads+unmapped_reads]:
        print('ERROR: check the star_logFinal_parser')
        print('if not in_read in [mapped_reads+unmapped_reads]:')
        print('in_read', in_read)
        print('mapped_reads', mapped_reads)
        print('mapped_reads_uniquely', mapped_reads_uniquely)
        print('mapped_reads_multiple', mapped_reads_multiple)
        print('unmapped_reads', unmapped_reads)
        sys.exit()

    stats_dic = dict()
    stats_dic.setdefault('mapped reads', mapped_reads)
    stats_dic.setdefault('unmapped reads', unmapped_reads)
    stats_dic.setdefault('uniquely mapped reads', mapped_reads_uniquely)
    return stats_dic

def star_SJout_parser(fn):
    spliced_reads = 0
    for line in open(fn):
        items = line.rstrip('\n').split('\t')
        spliced_uniquely_mapped_read = int(items[6])
        spliced_multi_mapped_read = int(items[7])
        spliced_reads += spliced_uniquely_mapped_read
        spliced_reads += spliced_multi_mapped_read

    stats_dic = dict()
    stats_dic.setdefault('spliced reads', spliced_reads)
    return stats_dic



def main(args):
    items = [args.sn]
    stat_dic = star_logFinal_parser(args.logfinal)
    items.append(stat_dic['mapped reads'])
    items.append(stat_dic['unmapped reads'])
    items.append(stat_dic['uniquely mapped reads'])
    stat_dic = star_SJout_parser(args.sjout)
    items.append(stat_dic['spliced reads'])

    with open(args.outfn, 'w') as outfh:
        items = [str(x) for x in items]
        outfh.write('{0}\n'.format('\t'.join(items)))

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('sn')
    parser.add_argument('outfn')
    parser.add_argument('logfinal')
    parser.add_argument('sjout')
    args = parser.parse_args()
    main(args)


