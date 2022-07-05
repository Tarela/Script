chrdict={}
checkfile = open("chr_limit_mm9.bed")
inf = open("/mnt/Storage/home/huse/ChongZhi/Data/peakfile/mouse_union_s150_c10.bed")
outf =open("/mnt/Storage/home/huse/ChongZhi/Data/peakfile/mouse_union_clap.bed",'w')
for line in checkfile:
    ll = line.split()
    if len(ll)==0:
	continue
    chrdict[ll[0]]=int(ll[2])
#print chrdict
for line in inf:
    ll = line.split()
    if int(ll[1])<1 or int(ll[2]) > chrdict[ll[0]]:
	print line
    else:
	outf.write(line)
checkfile.close()
inf.close()
outf.close()
