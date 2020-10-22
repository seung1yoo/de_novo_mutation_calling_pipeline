

import os
import sys
import json
import xlsxwriter

def words_to_outline_pe(words, infoDic):
    outline = list()
    TotalReadCnt = int(words[0])
    TotalLength  = int(words[1])
    TotalGCCnt   = int(words[2])
    NZeroReadCnt = int(words[3])
    N5ReadCnt    = int(words[4])
    TotalNCnt    = int(words[5])
    TotalQ30     = (int(words[8])  + int(words[9]))
    TotalQ20     = (int(words[12]) + int(words[13]))
    TotalLengthGB = "%.2f Gb" % (TotalLength * 1.0 / 1000000000)
    GCRate        = TotalGCCnt   * 100.0 / TotalLength
    NzRate        = NZeroReadCnt * 100.0 / TotalReadCnt
    N5Rate        = N5ReadCnt    * 100.0 / TotalReadCnt
    TotalNRate    = TotalNCnt    * 100.0 / TotalLength
    TotalQ30Rate  = TotalQ30     * 100.0 / TotalLength
    TotalQ20Rate  = TotalQ20     * 100.0 / TotalLength
    outline.append(infoDic['cst_id'])
    outline.append(infoDic['tbi_id'])
    outline.append(str(TotalReadCnt))
    outline.append(str(TotalLength))
    outline.append(str(TotalLengthGB))
    outline.append(str(TotalGCCnt))
    outline.append("%.2f%%" % GCRate)
    outline.append(str(NZeroReadCnt))
    outline.append("%.2f%%" % NzRate)
    outline.append(str(N5ReadCnt))
    outline.append("%.2f%%" % N5Rate)
    outline.append(str(TotalNCnt))
    outline.append("%.2f%%" % TotalNRate)
    outline.append(str(TotalQ30))
    outline.append("%.2f%%" % TotalQ30Rate)
    outline.append(str(TotalQ20))
    outline.append("%.2f%%" % TotalQ20Rate)
    return outline

def words_to_outline_se(words, infoDic):
    outline = list()
    TotalReadCnt = int(words[0])
    TotalLength  = int(words[1])
    TotalGCCnt   = int(words[2])
    NZeroReadCnt = int(words[3])
    N5ReadCnt    = int(words[4])
    TotalNCnt    = int(words[5])
    TotalQ30     = int(words[7])
    TotalQ20     = int(words[9])
    TotalLengthGB = "%.2f Gb" % (TotalLength * 1.0 / 1000000000)
    GCRate        = TotalGCCnt   * 100.0 / TotalLength
    NzRate        = NZeroReadCnt * 100.0 / TotalReadCnt
    N5Rate        = N5ReadCnt    * 100.0 / TotalReadCnt
    TotalNRate    = TotalNCnt    * 100.0 / TotalLength
    TotalQ30Rate  = TotalQ30     * 100.0 / TotalLength
    TotalQ20Rate  = TotalQ20     * 100.0 / TotalLength
    outline.append(infoDic['cst_id'])
    outline.append(infoDic['tbi_id'])
    outline.append(str(TotalReadCnt))
    outline.append(str(TotalLength))
    outline.append(str(TotalLengthGB))
    outline.append(str(TotalGCCnt))
    outline.append("%.2f%%" % GCRate)
    outline.append(str(NZeroReadCnt))
    outline.append("%.2f%%" % NzRate)
    outline.append(str(N5ReadCnt))
    outline.append("%.2f%%" % N5Rate)
    outline.append(str(TotalNCnt))
    outline.append("%.2f%%" % TotalNRate)
    outline.append(str(TotalQ30))
    outline.append("%.2f%%" % TotalQ30Rate)
    outline.append(str(TotalQ20))
    outline.append("%.2f%%" % TotalQ20Rate)
    return outline

def rst_to_xls(sample_dic, outfn, header): # modified other script
    outfh = open(outfn, 'w')
    if header in ['T']:
        headers = ['SampleID','TBI_ID','TotalReads','TotalBases','TotalBases(Gb)',
                   'GC_Count','GC_Rate','N_ZeroReads','N_ZeroReadsRate','N5_LessReads',
                   'N5_LessReadsRate','N_Count','N_Rate','Q30_MoreBases','Q30_MoreBasesRate',
                   'Q20_MoreBases','Q20_MoreBasesRate']
        outfh.write('{0}\n'.format('\t'.join(headers)))
    else:
        pass
    for statFile, infoDic in sorted(sample_dic.items()):
        for stats in open(statFile):
            if stats.startswith('#') :
                continue
            #
            words = stats.strip('\n').split("\t")
            if len(words) in [16]:
                outline = words_to_outline_pe(words, infoDic) 
            elif len(words) in [10]:
                outline = words_to_outline_se(words, infoDic) 
            else:
                print('ERROR : un-expected length of columns')
                sys.exit()
        outfh.write('{0}\n'.format('\t'.join(outline)))
    outfh.close()

def rst_to_xlsx(sample_dic, xls_fn, header): # modified other script
    if xls_fn.endswith('.xls'):
        xlsx_fn = '{0}x'.format(xls_fn)
    else:
        xlsx_fn = '{0}.xlsx'.format(xls_fn)

    book = xlsxwriter.Workbook(xlsx_fn)
    sheet1 = book.add_worksheet('Sequencing Result')
    colalign = book.add_format({'align' : 'center', 'valign' : 'center'})
    sheet1.set_column(0,0,13,   colalign)
    sheet1.set_column(1,1,13,   colalign)
    sheet1.set_column(2,2,11,   colalign)
    sheet1.set_column(3,3,13,   colalign)
    sheet1.set_column(4,4,15,   colalign)
    sheet1.set_column(5,5,11,   colalign)
    sheet1.set_column(6,6,9,    colalign)
    sheet1.set_column(7,7,13,   colalign)
    sheet1.set_column(8,8,18,   colalign)
    sheet1.set_column(9,9,14,   colalign)
    sheet1.set_column(10,10,19, colalign)
    sheet1.set_column(11,11,9,  colalign)
    sheet1.set_column(12,12,8,  colalign)
    sheet1.set_column(13,13,17, colalign)
    sheet1.set_column(14,14,20, colalign)
    sheet1.set_column(15,15,17, colalign)
    sheet1.set_column(16,16,20, colalign)

    sheet1header = [{'header' : 'SampleID'},
                    {'header' : 'TBI_ID'},
                    {'header' : 'TotalReads'},
                    {'header' : 'TotalBases'},
                    {'header' : 'TotalBases(Gb)'},
                    {'header' : 'GC_Count'},
                    {'header' : 'GC_Rate'},
                    {'header' : 'N_ZeroReads'},
                    {'header' : 'N_ZeroReadsRate'},
                    {'header' : 'N5_LessReads'},
                    {'header' : 'N5_LessReadsRate'},
                    {'header' : 'N_Count'},
                    {'header' : 'N_Rate'},
                    {'header' : 'Q30_MoreBases'},
                    {'header' : 'Q30_MoreBasesRate'},
                    {'header' : 'Q20_MoreBases'},
                    {'header' : 'Q20_MoreBasesRate'}]
    sheet1table = list()
    for statFile, infoDic in sorted(sample_dic.items()):
        for stats in open(statFile):
            if stats.startswith('#') :
                continue
            #
            words = stats.strip('\n').split("\t")
            if len(words) in [16]:
                outline = words_to_outline_pe(words, infoDic) 
            elif len(words) in [10]:
                outline = words_to_outline_se(words, infoDic) 
            else:
                print('ERROR : un-expected length of columns')
                sys.exit()
            outline = [str(x) for x in outline]
            sheet1table.append(outline)

    if header in ['T']:
        sheet1.add_table(0,0,(len(sample_dic)),16,{'data' : sheet1table,
                                            'columns' : sheet1header,
                                            'style': 'Table Style Light 15',
                                            'autofilter': False,
                                            'banded_rows': 0})
    else:
        sheet1.add_table(0,0,(len(sample_dic)-1),16,{'data' : sheet1table,
                                            'header_row' : False,
                                            'style': 'Table Style Light 15',
                                            'autofilter': False,
                                            'banded_rows': 0})

    book.close()


def main(args):

    rst_to_xls(json.load(open(args.json)), args.outfn, args.header)
    rst_to_xlsx(json.load(open(args.json)), args.outfn, args.header)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('json')
    parser.add_argument('outfn')
    parser.add_argument('--header', choices=('T','F'), default='T')
    args = parser.parse_args()
    main(args)
