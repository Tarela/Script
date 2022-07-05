#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description:
Script for plot paired end size distribution


"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys
from optparse import OptionParser

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def fraglen(inputfile,output,toplim,bottomlim):
  #  output += 'out_'
    inf = open(inputfile)
    D = {}
    #print int(float(bottomlim)),int(float(toplim))
    for i in range(int(float(bottomlim)),int(float(toplim))):
        D[i]=0
    
    for line in inf:
        if line.startswith('@'):
            continue
        ll = line.strip().split("\t")
        if ll[2].startswith('chr') and int(ll[8]) > 0 :
            if D.has_key(int(ll[8])):
                D[int(ll[8])]+=1    
    inf.close()
    length = []
    times = []
    for i in sorted(D.keys()):
        length.append(i)
        times.append(D[i])
    outf = open(output+'.r','w')
    outf.write('fraglen<-c('+str(length)[1:-1]+')\n')
    outf.write('%s<-c('%(output)+str(times)[1:-1]+')\n')
    plotscript = """
pdf(file="%s.pdf")
plot(fraglen,%s,type="l",lwd=2,col="blue",xlab="Pairend size",ylab="count",main="Paired end size distribution of %s")
legend("topleft",c("%s"),col=c("blue"),lwd=5,bty="n")
dev.off()
"""%(output,output,output,output)
    outf.write(plotscript)
    outf.close()
    os.system('Rscript %s.r'%(output))


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog -i inputfile -o outname [options] "
    description = """Script for plot paired end size distribution , input file should be bowtie paired end mapping result in .sam format"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="input file in paired end sam format , recommand is bowtie1 paired end mapping result")
    optparser.add_option("-o","--outname",dest="output",type="str",default='NA',
                         help="name of your output file, you will have 2 output file named as youname.r and yourname.pdf")

#========minor options=============
    optparser.add_option("--bottomlimit",dest="botlim",type="int",default =0,
                         help="bottom limit of your paired end size range, default is 0")

    optparser.add_option("--upperlimit",dest="uplim",type="int",default =600,
                         help="upper limit of your paired end size range, default is 600")

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    
    if not inputfile or not output:
        optparser.print_help()
        sys.exit(1)
    
    fraglen(inputfile,output,options.uplim,options.botlim)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

