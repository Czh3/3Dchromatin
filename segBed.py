import sys
window = 1000
for i in open(sys.argv[1]):
	j = i.strip().split()
	for k in range(int(j[1]), int(j[2]), window):
		print j[0]+"\t"+str(k)+"\t"+str(k+window)
