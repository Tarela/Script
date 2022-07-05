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
    FLdict = [0]*uselim
    for line in inf:
        ll = line.strip().split("\t")
        fraglen = int(ll[2]) - int(ll[1])
        if fraglen > uselim or fraglen <= 0:
            print inputfile, ll
            continue
        FLdict[(fraglen-1)] += 1
    inf.close()
    outf = open(output,'w')
    outf.write("\t".join(map(str,FLdict))+"\n")
    outf.close()


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
                         help="name of your output file")

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

