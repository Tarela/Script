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
ymax<-%s
xmax<-%s
pdf(file = "%s")
plot(plus,type="l",col="darkred",xlim=c(1,xmax),ylim=c(0,ymax))
par(new=T)
plot(minus,type="l",col="darkblue",xlim=c(1,xmax),ylim=c(0,ymax))
legend('topleft',c('plus','minus'),col=c('darkred','darkblue'),lwd=3)
dev.off()
'''%(str(plus)[1:-1],str(minus)[1:-1],str(max(plus+minus)),len(plus),out+'.pdf')
    Rscript.write(Rstring)
    Rscript.close()
    os.system('Rscript %s'%(out+'.r'))

def make_template(data,flank,topmotif,out):
    i =0
    datalist = []
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
        datalist.append(ll)
        if ll[814] == '1' :
            templatelist.append(ll)
    inf.close()
    templatelist.sort(key = lambda x:float(x[4]),reverse=True)
    for ll in templatelist[:topmotif]:
        if ll[5] == "+":
            pp.append(ll[14:414][(201 - flank):(201 + ml + flank)])
            pm.append(ll[414:814][(201 - flank):(201 + ml + flank)])
        if ll[5] == '-' :
            pm.append(ll[14:414][::-1][(201 - ml - flank) : (201 -ml +ml +flank)])
            pp.append(ll[414:814][::-1][(201 - ml - flank) : (201 -ml +ml +flank)])

    meanp = apply_mean(pp)
    meanm = apply_mean(pm)
    allsum = sum(meanp)+sum(meanm)
    P=[]
    M=[]
    for i in range(len(meanp)):
        P.append(meanp[i]/allsum)
        M.append(meanm[i]/allsum)    

    plot_template(P,M,out)

    return (P,M,datalist,ml)    

def match_pattern(pattern_plus,pattern_minus,observe_plus,observe_minus,l,persudo):
    ### l is the length of pattern
    observe_sum = sum(observe_plus)+sum(observe_minus)
    if observe_sum == 0 :
        return "no"
    p0_sum = observe_sum + persudo*l
    lambda_plus = pattern_plus
    lambda_minus =  pattern_minus
    lambda_plus0 = [(x*1.0+persudo)/p0_sum for x in observe_plus]
    lambda_minus0 = [(x*1.0+persudo)/p0_sum for x in observe_minus]
    f1=0
    for i in range(l):
        f1 = f1 + observe_plus[i]*(math.log10(lambda_plus[i])-math.log10(lambda_plus0[i])) + observe_minus[i]*(math.log10(lambda_minus[i])-math.log10(lambda_minus0[i]))
    return (f1/observe_sum)
def scan(plus,minus,datalist,flank,out,ml):
    l = ml + 2*flank
    ii =0
    outf = open(out + '.txt','w')
    for ll in datalist:
        if ll[5] == "+":
            pob = (ll[14:414][(201 - flank):(201 + ml + flank)])
            mob = (ll[414:814][(201 - flank):(201 + ml + flank)])
        if ll[5] == '-' :
            mob = (ll[14:414][::-1][(201 - ml - flank) : (201 -ml +ml +flank)])
            pob = (ll[414:814][::-1][(201 - ml - flank) : (201 -ml +ml +flank)])
        #score = match_pattern(plus,minus,map(float,pob),map(float,mob),l,1)
        newll = ll[:6]+[ll[814],float(ll[6])+float(ll[7]),float(ll[8])+float(ll[9]),float(ll[10])+float(ll[11]),float(ll[12])+float(ll[13])]
        for i in range(1,21):
            newll.append(match_pattern(plus,minus,map(float,pob),map(float,mob),l,i*1.0/10))
        outf.write("\t".join(map(str,newll))+"\n")

    outf.close()
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
#========minor options=============
    optparser.add_option("--flank",dest="flank",type="int",default = 3,
                         help="+- N for flanking region , length of your template will be motiflength + 2N , default N = 3")
    optparser.add_option("--topmotif",dest="topmotif",type="int",default = 5000,
                         help="motif number for making tamplate, default is 5000")


    (options,args) = optparser.parse_args()

    data = options.data
    out = options.out
    flank = options.flank
    topmotif = options.topmotif

    if not data :
        optparser.print_help()
        sys.exit(1)
    
    (P,M,dalist,ml) = make_template(data,flank,topmotif,out)
    scan(P,M,dalist,flank,out,ml)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

