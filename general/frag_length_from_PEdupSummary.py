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

def fraglen(inputfile,output,LIM):
  #  output += 'out_'
    uselim = int(float(LIM))
    inf = open(inputfile)
    DnoMall = []
    DnoMuniq = []
    DchrMall = []
    DchrMuniq = []
    #print int(float(bottomlim)),int(float(toplim))
    for i in range(uselim):
        DnoMall.append(0)
        DnoMuniq.append(0)
        DchrMall.append(0)
        DchrMuniq.append(0)
    
    for line in inf:
        ll = line.strip().split("\t")
        fraglen = int(ll[8])
        dupN = int(ll[4])
        if fraglen > uselim or fraglen == 0:
            continue
        if ll[0] == "chrM":
            DchrMuniq[(fraglen-1)] += 1
            DchrMall[(fraglen-1)] += dupN
        else:
            DnoMuniq[(fraglen-1)] += 1
            DnoMall[(fraglen-1)] += dupN
    inf.close()
    outf = open(output+'_fraglen.txt','w')
    newll1 = [output+"_noMall"] + DnoMall
    newll2 = [output+"_noMuniq"] + DnoMuniq
    newll3 = [output+"_chrMall"] + DchrMall
    newll4 = [output+"_chrMuniq"] + DchrMuniq
    outf.write("\t".join(map(str,newll1))+"\n")
    outf.write("\t".join(map(str,newll2))+"\n")
    outf.write("\t".join(map(str,newll3))+"\n")
    outf.write("\t".join(map(str,newll4))+"\n")
    


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

    optparser.add_option("--limit",dest="limit",type="int",default =2000,
                         help="upper limit of your paired end size range, default is 2000")

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    
    if not inputfile or not output:
        optparser.print_help()
        sys.exit(1)
    
    fraglen(inputfile,output,options.limit)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

