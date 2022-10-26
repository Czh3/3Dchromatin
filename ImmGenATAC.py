import sys


bed = {}
for i in open(sys.argv[1]):
	j = i.strip().split("\t")
	pos1 = int(j[1])
	pos2 = int(j[2])
	if bed.has_key(j[0]):
		bed[j[0]].append([pos1, pos2])
	else:
		bed[j[0]] = [[pos1, pos2]]

OCR = {}
for i in open("/lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/cell_2019/ImmGenATAC18_AllOCRsInfo.tab").readlines()[1:]:
	i = i.strip()
	j = i.split("\t")
	chrom = j[1]
	summit = int(j[2])
	if OCR.has_key(chrom):
		OCR[chrom].append([summit, j[8:]])
	else:
		OCR[chrom] = [[summit, j[8:]]]

# output
for i in bed:
	for j in OCR:
		if i == j:
			for ii in bed[i]:
				for jj in OCR[j]:
					if jj[0] in xrange(ii[0], ii[1]):
						print "\t".join(jj[1:])


