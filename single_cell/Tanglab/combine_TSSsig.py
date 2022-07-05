import os,sys
def addfile(usedict,add_file):
    inf = open(add_file)
    for line in inf:
        ll = line.split()
        name = ll[3]
        ov = ll[-1]
        if usedict.has_key(name):
            usedict[name].append(ov)
        else:
            usedict[name] = [ov]



TSS = {}

CN = []
ATAConly = "/scratch/sh8tv/Project/scATAC/Data/scATAC/Tanglab_ATAC/process/ATAConly/"
coassay = "/scratch/sh8tv/Project/scATAC/Data/scATAC/Tanglab_ATAC/process/coassay_ATAC/"
allf = sorted(os.listdir(ATAConly))
for f in allf:
    if f.startswith("TSS3kb_m") and f.endswith("allreads.bed"):
        thisname = f[7:-13]
        CN.append(thisname)
        TSSfile = ATAConly+f
#        enhancerfile = "ESCenhancer_%s_allreads.bed"%(thisname)
        addfile(TSS,TSSfile)
#        addfile(enhancer,enhancerfile)
allf = sorted(os.listdir(coassay))
for f in allf:
    if f.startswith("TSS3kb_m") and f.endswith("allreads.bed"):
        thisname = f[7:-13]
        CN.append(thisname)
        TSSfile = coassay+f
#        enhancerfile = "ESCenhancer_%s_allreads.bed"%(thisname)
        addfile(TSS,TSSfile)




outfTSS = open("TSS3kb_allreads_overlap.txt",'w')
#outfenhancer = open("ESCenhancer_allreads_overlap.txt",'w')

outfTSS.write("\t".join(CN)+"\n")
#outfenhancer.write("\t".join(CN)+"\n")

for gname in sorted(TSS.keys()):
    outfTSS.write("\t".join([gname]+TSS[gname])+"\n")

#for pname in sorted(enhancer.keys()):
#    outfenhancer.write("\t".join([pname]+enhancer[pname])+"\n")

outfTSS.close()
#outfenhancer.close()