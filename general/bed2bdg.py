#!/usr/bin/env python
import sys
inf = open(sys.argv[1])
outf = open(sys.argv[2],'w')
chrinfo ={}
Chrinfo = open('/nv/vol190/zanglab/sh8tv/Data/Genome/hg38/hg38_clean.len')
for line in Chrinfo:
    ll = line.split()
    chrinfo[ll[0]] = int(ll[1])
Chrinfo.close()

nowchr = 'na'
for line in inf: 
    ll = line.split()
    if not ll[0] in chrinfo.keys():
        continue
    if nowchr != ll[0]:
        if nowchr == 'na':
            pass
        else:
            if not nowloci == chrinfo[nowchr]:
                pass#newll = [nowchr,nowloci,chrinfo[nowchr],-1]
                #outf.write("\t".join(map(str,newll))+"\n")
            
        nowchr = ll[0]
        #newll = [nowchr,0,ll[1],-1]
        #outf.write("\t".join(map(str,newll))+"\n")
        newll = [nowchr,ll[1],int(ll[1])+1,ll[3]]
        outf.write("\t".join(map(str,newll))+"\n")
        nowloci = int(ll[1])+1
    else:
        #newll = [nowchr,nowloci,ll[1],-1]
        #outf.write("\t".join(map(str,newll))+"\n")
        newll = [nowchr,ll[1],int(ll[1])+1,ll[3]]
        outf.write("\t".join(map(str,newll))+"\n")
        nowloci = int(ll[1])+1

if not nowloci == chrinfo[nowchr]:
    pass#newll = [nowchr,nowloci,chrinfo[nowchr],-1]
    #outf.write("\t".join(map(str,newll))+"\n")

inf.close()
outf.close()
