
import math

class DEGParser:
    def __init__(self, fn, cn):
        self.comp_name = cn
        self.ctrl_id = ''
        self.case_id = ''
        self.deg_dic = self.make_dic(fn)

    def make_dic(self, fn):
        deg_dic = dict()
        for line in open(fn):
            items = line.rstrip('\n').split('\t')
            #['test_id', 'gene_id', 'gene', 'locus', 'sample_1', 'sample_2', 'status', 'value_1', 'value_2', 'log2(fold_change)', 'test_stat', 'p_value', 'q_value', 'significant']
            if items[0] in ['test_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            #['0610005C13Rik', '0610005C13Rik', '0610005C13Rik', 'chr7:45567794-45589710', 'AKTPF_LOH_Acv', 'AKTPF_hez_Acv', 'NOTEST', '0', '0.0309635', 'inf', '0', '1', '1', 'no']
            gid = items[idx_dic['test_id']]
            self.ctrl_id = items[idx_dic['sample_1']]
            self.case_id = items[idx_dic['sample_2']]
            status = items[idx_dic['status']]
            ctrl_value = items[idx_dic['value_1']]
            case_value = items[idx_dic['value_2']]
            log2fc = math.log(float(case_value)+1.0,2) - math.log(float(ctrl_value)+1.0,2)
            p = items[idx_dic['p_value']]
            q = items[idx_dic['q_value']]

            deg_dic.setdefault(gid, {}).setdefault('status', status)
            deg_dic.setdefault(gid, {}).setdefault('ctrl_id', self.ctrl_id)
            deg_dic.setdefault(gid, {}).setdefault('case_id', self.case_id)
            deg_dic.setdefault(gid, {}).setdefault('ctrl_value', ctrl_value)
            deg_dic.setdefault(gid, {}).setdefault('case_value', case_value)
            deg_dic.setdefault(gid, {}).setdefault('log2fc', log2fc)
            if log2fc > 0.0:
                deg_dic.setdefault(gid, {}).setdefault('updown', 'up')
            elif log2fc < 0.0:
                deg_dic.setdefault(gid, {}).setdefault('updown', 'down')
            else:
                deg_dic.setdefault(gid, {}).setdefault('updown', 'flat')
            deg_dic.setdefault(gid, {}).setdefault('p', p)
            deg_dic.setdefault(gid, {}).setdefault('q', q)
            if abs(log2fc) >= 1 and float(p) < 0.05:
                deg_dic.setdefault(gid, {}).setdefault('isdeg', 'Y')
            else:
                deg_dic.setdefault(gid, {}).setdefault('isdeg', 'N')


        return deg_dic

def main():
    deg01 = DEGParser('/data06/project/VIP_RNASeq_20200529/analysis/cuffdiff/DEG01/gene_exp.diff', 'DEG01')
    deg02 = DEGParser('/data06/project/VIP_RNASeq_20200529/analysis/cuffdiff/DEG02/gene_exp.diff', 'DEG02')
    deg03 = DEGParser('/data06/project/VIP_RNASeq_20200529/analysis/cuffdiff/DEG03/gene_exp.diff', 'DEG03')
    deg04 = DEGParser('/data06/project/VIP_RNASeq_20200529/analysis/cuffdiff/DEG04/gene_exp.diff', 'DEG04')

    comp_s = ['DEG01','DEG02','DEG03','DEG04']
    obj_s = [deg01, deg02, deg03, deg04]
    g_dic = dict()
    for obj in obj_s:
        print(obj.comp_name, obj.ctrl_id, obj.case_id)
        for gid, info_dic in obj.deg_dic.items():
            g_dic.setdefault(gid, {}).setdefault(obj.comp_name, info_dic)

    outfh = open('/data06/project/VIP_RNASeq_20200529/analysis/cuffdiff/cuffdiff.xls','w')
    headers = ['gene_id']
    for comp_name in comp_s:
        for prefix in ['status','ctrl_id','case_id','ctrl_value','case_value','log2fc','p','q','updown','isdeg']:
            header = '{0}:{1}'.format(comp_name, prefix)
            headers.append(header)
    outfh.write('{0}\n'.format('\t'.join(headers)))
    for g_id, comp_dic in g_dic.items():
        items = [g_id]
        for comp_name in comp_s:
            info_dic = comp_dic[comp_name]
            items.append(info_dic['status'])
            items.append(info_dic['ctrl_id'])
            items.append(info_dic['case_id'])
            items.append(info_dic['ctrl_value'])
            items.append(info_dic['case_value'])
            items.append(str(round(info_dic['log2fc'], 2)))
            items.append(info_dic['p'])
            items.append(info_dic['q'])
            items.append(info_dic['updown'])
            items.append(info_dic['isdeg'])
        outfh.write('{0}\n'.format('\t'.join(items)))
    outfh.close()





if __name__=='__main__':
    main()
