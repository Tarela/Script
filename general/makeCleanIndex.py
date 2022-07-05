chrmlist = []
inf = open('chromInfo_mm10_clean.txt')
for line in inf: 
    chrmlist.append(line.split()[0])

inf = open('mm10.fa')
outf = open('mm10_clean.fa','w')
for line in inf:
    if line.startswith(">"):
        if line.strip()[1:] in chrmlist:
            WRITE=1
        else:
            WRITE=0
    if WRITE == 1:
        outf.write(line)
outf.close()
inf.close()

