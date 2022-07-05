import sys
infile = sys.argv[1]
eachcount =  int(sys.argv[2])
extention = sys.argv[3]
inf = open(infile)
each = eachcount
b=0
c=0
outf = open("tmp_%s.%s"%(str(c/each),extention),'w')
for line in inf:
    outf.write(line)
    c+=1
    if c%each==0:
        outf.close()
        outf=open("tmp_%s.%s"%(str(c/each),extention),'w')

outf.close()
