'''
Created on 2011-5-24

@author: undergraduate
'''
#!/usr/bin/env python
#Time-stamp:<2011-05-24 Sheng'en>
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

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def bed2fa(bedfile,assembly):
    cmd = "fastaFromBed -fi %s -bed %s -fo %s -name"%(assembly,bedfile,"/tmp/tmp123.tmp")
    os.system(cmd)

def CG(bedfile,assembly,out):
    bed2fa(bedfile,assembly)
    beddict={}
    cgdict={}
    inf = open(bedfile)
    for line in inf:
	ll = line.split()
	beddict[ll[3]]=ll
    inf.close()
    inf = open("/tmp/tmp123.tmp")
    for line in   inf:
	if line.startswith(">"):
	    name = line.strip().lstrip(">")
	else:
	    seq = line.strip()
	    cg=0
	    length=0
	    for bp in seq:
		length+=1
		if bp=="C"or bp=="G":
		    cg+=1
	    ratio=str(cg*1.0/length)
	    cgdict[name]=[ratio]
    outf = open(out,'w')
    for name in beddict:
	newll=beddict[name]+cgdict[name]
	newline = "\t".join(newll)+"\n"
	outf.write(newline)
    outf.close()
    os.system("rm %s"%("/tmp/tmp123.tmp"))



# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog <-b bedfile> <-a assembly> <-o out> "
    description = """cg content for each region"""

    optparser = OptionParser(version="%prog 0.72",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========core options=============
    optparser.add_option("-b","--bedfile",dest="bedfile",type="str",
                         help="")
    optparser.add_option("-a","--assembly",dest="assembly",type="str",
                         help="")
    optparser.add_option("-o","--outfile",dest="outfile",type="str",
                         help="")
    (options,args) = optparser.parse_args()

    bedfile = options.bedfile
    assembly = options.assembly
    outfile = options.outfile
    if not bedfile:
        optparser.print_help()

        sys.exit(1)
    
    CG(bedfile,assembly,outfile)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

