'''
Created on XXXX-XX-XX

@author: Tarela
'''
#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging
import string
import math
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def apply_mean(vector):
    mv = [0]*len(vector[0])
    for v in vector:
        for i in range(len(mv)):
            mv[i] += float(v[i])
    for i in range(len(mv)):
        mv[i] = mv[i]*1.0/len(vector)

    return mv
def plot_template(plus,minus,out):
    Rscript = open(out+'.r','w')
    Rstring = '''
plus<-c(%s)
minus<-c(%s)
xmax<-%s
pdf(file = "%s")
plot(plus,type="l",col="darkred",xlim=c(0,xmax),lwd=2,xlab="Relative distance (bp)",ylab="Cleavage(left),Conserve(right)",main="",axes=F)
axis(side=2)
par(new=T)
plot(minus,type="l",col="darkblue",xlim=c(0,xmax),lwd=2,xlab="",ylab="",main="",axes=F)
axis(side=4)
abline(v=xmax/2,lwd=2,lty=2)
legend('topleft',c('Cleavage','Conserve'),col=c('darkred','darkblue'),lwd=3,bty="n")
axis(side=1,at = seq(0,xmax,10),labels=seq(-xmax/2,xmax/2,10))
box()
dev.off()
'''%(str(plus)[1:-1],str(minus)[1:-1],len(plus),out+'.pdf')
    Rscript.write(Rstring)
    Rscript.close()
    os.system('Rscript %s'%(out+'.r'))

def make_template(data,flank,pflank,topmotif,out,pbw,mbw):
    w_plus_H=BigWigFile(open(pbw, 'rb'))
    w_minus_H=BigWigFile(open(mbw, 'rb'))
    i =0
    templatelist = []
    pp=[]
    pm=[]
    inf = open(data)
    l1st = inf.readline().split()
    ml = int(l1st[2])-int(l1st[1])
    inf.seek(0)
    for line in inf:
        #if i >= topmotif:
         #   break
        ll = line.split()
        templatelist.append(ll)

    inf.close()
    templatelist.sort(key = lambda x:float(x[4]),reverse=True)
    for ll in templatelist:
        p_sum = list(w_plus_H.summarize(ll[0],int(ll[1])-flank,int(ll[1])+flank,2*flank).sum_data)
        m_sum = list(w_minus_H.summarize(ll[0],int(ll[1])-flank,int(ll[1])+flank,2*flank).sum_data)
        if 1:
            pp.append(p_sum[(flank + 1 + ml/2 -pflank):(flank +1 + ml/2  + pflank)])
            pm.append(m_sum[(flank +1 + ml/2 - pflank):(flank +1 + ml/2  + pflank)])
        #if ll[5] == '-' :
        #    pm.append(p_sum[::-1][(flank +1 + ml/2 - ml - pflank) : (flank +1 + ml/2 -ml +pflank)])
        #    pp.append(m_sum[::-1][(flank +1 + ml/2 - ml - pflank) : (flank +1 + ml/2 -ml +pflank)])

    meanp = apply_mean(pp)
    meanm = apply_mean(pm)
    allsum = sum(meanp)+sum(meanm)
    P=[]
    M=[]
    for i in range(len(meanp)):
        P.append(meanp[i])#/allsum)
        M.append(meanm[i])#/allsum)    

    plot_template(P,M,out)

   
# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputdata",dest="data",type="str",
                         help="")
    optparser.add_option("-o","--outputdata",dest="out",type="str",
                         help="")
    optparser.add_option("-p","--pbw",dest="pbw",type="str",
                         help="")
    optparser.add_option("-m","--mbw",dest="mbw",type="str",
                         help="")
#========minor options=============
    optparser.add_option("--flank",dest="flank",type="int",default = 100,
                         help="+- N for flanking region , default +- 100")
    optparser.add_option("--pflank",dest="pflank",type="int",default = 50,
                         help="+- N for plot flanking region default +- 50")


    optparser.add_option("--topmotif",dest="topmotif",type="int",default = 5000,
                         help="motif number for making tamplate, default is 5000")


    (options,args) = optparser.parse_args()

    peak = options.data
    out = options.out
    flank = options.flank
    topmotif = options.topmotif
    pflank = options.pflank
    pbw = options.pbw
    mbw = options.mbw
    

    if not peak :
        optparser.print_help()
        sys.exit(1)
    
    make_template(peak,flank,pflank,topmotif,out,pbw,mbw)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


