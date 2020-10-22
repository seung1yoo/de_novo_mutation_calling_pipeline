#!/usr/bin/python
import sys, gzip, subprocess, getopt, os, commands

PLATFORM_ZERO_DIC = {"Sanger" : -33, "Solexa" : -64, "Illumina 1.3+" : -64, "Illumina 1.5+" : -64, "Illumina 1.8+" : -33}

FIRST_IN_FASTQ = None
RESULT_OUTPUT = None

def init():
	global FIRST_IN_FASTQ, SECOND_IN_FASTQ, RESULT_OUTPUT
	options, args = getopt.getopt(sys.argv[1:],"")
	
	FIRST_IN_FASTQ = args[0]
	RESULT_OUTPUT = args[1]

def checkFlatform():
        # quality #
        qualityIdcDic = {}
        firReadFile = gzip.open(FIRST_IN_FASTQ, "rb")
        lineIdx = -1
        for line in firReadFile :
                lineIdx += 1
                if lineIdx == 100000 : break
                if lineIdx % 4 != 3 : continue
                for idc in line.strip() : qualityIdcDic[idc] = None
        firReadFile.close()

        platform = None
        qualityIdcDicKeys = qualityIdcDic.keys()
        minQualityIdc = min(qualityIdcDicKeys)
        maxQualityIdc = max(qualityIdcDicKeys)
        if   minQualityIdc == "!" and maxQualityIdc == "I" :
                platform = "Sanger"
        elif minQualityIdc == ";" and maxQualityIdc == "h" :
                platform = "Solexa"
        elif minQualityIdc == "@" and maxQualityIdc == "h" :
                platform = "Illumina 1.3+"
        elif minQualityIdc == 'B' and maxQualityIdc == "h" :
                platform = "Illumina 1.5+"
        elif minQualityIdc == "!" and maxQualityIdc == "J" :
                platform = "Illumina 1.8+"
        else :
                if   minQualityIdc < ";" :
                        platform = "Illumina 1.8+"
                elif minQualityIdc < "@" :
                        platform = "Solexa"
                elif minQualityIdc < "B" :
                        platform = "Illumina 1.3+"
                else :
                        platform = "Illumina 1.5+"
        return platform

# N and CG count
def NandGCCount(_seq) :
	N_count = 0
	GC_count = 0
	countLen = len(_seq)
	for s in _seq :
		if s == "N" :
			N_count += 1
		elif s == "G" or s == "C" :
			GC_count += 1
		else :
			pass
	return N_count, GC_count

# qulity count
def QualityCount(_qualities, _flatform) :
	convertValue = PLATFORM_ZERO_DIC[_flatform]
	qualitiesLen = len(_qualities)
	qualityQ20MoreTotal = 0
	qualityQ30MoreTotal = 0
	qualityQ20MoreRead = 0
	qualityQ30MoreRead = 0
	TotalQValue = 0
	for quality in _qualities :
		qValue = ord(quality) + convertValue
		TotalQValue += qValue
		if qValue >= 20 :
			qualityQ20MoreTotal += 1
		if qValue >= 30 :
			qualityQ30MoreTotal += 1
	if 30 <= (TotalQValue * 1.0 / qualitiesLen) :
		qualityQ30MoreRead += 1
	if 20 <= (TotalQValue * 1.0 / qualitiesLen) :
		qualityQ20MoreRead += 1
	return qualityQ20MoreTotal, qualityQ30MoreTotal, qualityQ20MoreRead, qualityQ30MoreRead

if __name__ == "__main__" :
	if len(sys.argv) != 3 :
		print "Usage : python %s <in.first.fq.gz> <out.result.txt>" % sys.argv[0]
		sys.exit()
	else :
		init()

	flatform = checkFlatform()

	fiFirst = gzip.open(FIRST_IN_FASTQ, "rb")
#	fiFirst = os.popen("zcat %s" % FIRST_IN_FASTQ)

	totalReadCnt = 0
	totalLen = 0
	totalGCCnt = 0
	totalQ30AboveReadCnt = 0
	totalQ20AboveReadCnt = 0
	totalQ30AboveBaseCnt = 0
	totalQ20AboveBaseCnt = 0
	totalNCnt = 0
	totalNZeroCnt = 0
	totalN5BelowCnt = 0

	while True :
		line1_id = fiFirst.readline().rstrip()
		line1_seq = fiFirst.readline().rstrip()
		line1_strand = fiFirst.readline().rstrip()
		line1_quality = fiFirst.readline().rstrip()

		if not line1_id :
			break

		totalReadCnt += 1
		totalLen += len(line1_seq)
		nCount, gcCount = NandGCCount(line1_seq)
		if nCount == 0 :
			totalNZeroCnt += 1
		if (nCount * 1.0 / len(line1_seq)) < 0.05 :
			totalN5BelowCnt += 1
		totalNCnt += nCount
		totalGCCnt += gcCount
		
		check_QB20, check_QB30, check_QR20, check_QR30 = QualityCount(line1_quality, flatform)
		totalQ30AboveReadCnt += check_QR30
		totalQ20AboveReadCnt += check_QR20
		totalQ30AboveBaseCnt += check_QB30
		totalQ20AboveBaseCnt += check_QB20

	fo = open(RESULT_OUTPUT, "w")
	headBufList = []
	headBufList.append("totalReadCnt")
	headBufList.append("totalLength")
	headBufList.append("totalGCCnt")
	headBufList.append("NZeroIncludeReadCnt")
	headBufList.append("N5IncludeReadCnt")
	headBufList.append("totalNCnt")
	headBufList.append("totalQ30ReadCnt")
	headBufList.append("totalQ30BaseCnt")
	headBufList.append("totalQ20ReadCnt")
	headBufList.append("totalQ20BaseCnt")
	fo.write("#%s\n" % "\t".join(headBufList))

	valueBufList = []
	valueBufList.append(str(totalReadCnt))
	valueBufList.append(str(totalLen))
	valueBufList.append(str(totalGCCnt))
	valueBufList.append(str(totalNZeroCnt))
	valueBufList.append(str(totalN5BelowCnt))
	valueBufList.append(str(totalNCnt))
	valueBufList.append(str(totalQ30AboveReadCnt))
	valueBufList.append(str(totalQ30AboveBaseCnt))
	valueBufList.append(str(totalQ20AboveReadCnt))
	valueBufList.append(str(totalQ20AboveBaseCnt))
	fo.write("%s\n" % "\t".join(valueBufList))
	fo.close()

