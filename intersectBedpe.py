import glob
import os



def intersectBedpe(bedpe1, bedpe2):
	res = []
	for a in bedpe1:
		flag = 0
		for b in bedpe2:
			if a[0] == b[0] and (a[1] <= b[1] and  a[2] >= b[1]) and (a[4] <= b[4] and a[5] >= b[5]):
				res.append(a + [b[6]])
				flag = 1
				break
		if flag == 0:
			res.append(a + ["NA"])
	return(res)


bedpe = []
for i in open("../analysis/loops/merge.loops1.bedpe").readlines()[10:]:
	j = i.split()
	j[1],j[2],j[4],j[5] = int(j[1]), int(j[2]), int(j[4]), int(j[5])
	bedpe.append(j[0:6])

for i in glob.glob("/lustre/user/liclab/zhangc/proj/blood3D_3/analysis/loops/merge/*"):
	name = os.path.basename(i)
	os.system("cat %s/requested_list_25000.bedpe %s/requested_list_10000.bedpe %s/requested_list_5000.bedpe > %s/3requested_list.bedpe" % (i, i,i, i ))
	req = []
	for q in open(i+"/3requested_list.bedpe"):
		if q.startswith("chr1"):continue
		j = q.split()
		j[1],j[2],j[4],j[5] = int(j[1]), int(j[2]), int(j[4]), int(j[5])
		req.append(j[0:6] + [max([float(p) for p in j[20:]])])

	bedpe = intersectBedpe(bedpe, req)

for i in bedpe:
	print "\t".join([str(j) for j in i])



