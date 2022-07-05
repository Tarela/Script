#!/usr/bin/env python
import subprocess
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
rdep =  ["24cell","256cell","1Kcell","oblong","dome","shield","bud","28hpf","2dpf","5dpf"]
zygotic = ["total_128","total_256","total_512","zygotic_128_1","zygotic_128_2","zygotic_256_1","zygotic_256_2","zygotic_512_1","zygotic_512_2"]
expdict = {}
for sample in rdep:
    filename = '/mnt/Storage/home/wangwen/work2/OCR/public_data/RNA/10stages/'+sample+'/cuff_ref/isoforms.fpkm_tracking'
    print sp('wc -l %s'%filename)
    inf = open(filename)
    for line in inf:
        if line.startswith('tracking_id'):
            continue
        ll = line.strip().split("\t")
        gname = ll[0]
        loci = ll[6]
        exp = ll[9]
        KEY = gname+":"+loci
        if expdict.has_key(KEY):
            expdict[KEY].append(exp)
        else:
            expdict[KEY] = [exp]

for sample in zygotic:
    filename = '/mnt/Storage/home/wangwen/work2/OCR/public_data/RNA/GSE47709/cuff_ref/'+sample+'/isoforms.fpkm_tracking'
    print sp('wc -l %s'%filename)
    inf = open(filename)
    for line in inf:
        if line.startswith('tracking_id'):
            continue
        ll = line.strip().split("\t")
        gname = ll[0]
        loci = ll[6]
        exp = ll[9]
        KEY = gname+":"+loci
        if expdict.has_key(KEY):
            expdict[KEY].append(exp)
        else:
            expdict[KEY] = [exp]


symboldict = {}
inf = open('/mnt/Storage/home/huse/Data/refseq/zv9_refgene.txt')
for line in inf:
    ll = line.split()
    if "_" in ll[2]:
        continue
    if not symboldict.has_key(ll[1]):
        symboldict[ll[1]] = [ll[2],ll[4],ll[5],ll[3],ll[12]]
inf.close()


outf = open('zv9_10stages_exptable.txt','w')
newll = ['chrm','start','end','gname','symbol','strand','rowname'] + rdep + zygotic
outf.write("\t".join(newll)+"\n")
for i in sorted(expdict.keys()):
    gname = i.split(":")[0]
    if len(gname.split("_")) > 2:
        GN = "_".join(gname.split("_")[:2])
    else:
        GN = gname
    if symboldict.has_key(GN):
        symbol = symboldict[GN][4]
        chrom = symboldict[GN][0]
        start = symboldict[GN][1]
        end = symboldict[GN][2]
        strand = symboldict[GN][3]
        newll = [chrom,start,end,GN,symbol,strand,i] + expdict[i]
        outf.write("\t".join(newll)+"\n")
    else:
        symbol = '-'
        print GN
    
outf.close()

