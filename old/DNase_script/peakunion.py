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

def fetchpeak(table,directory,cutoff):
    if not directory.endswith("/"):
    	directory+="/"
    inf = open(table,'r')
#    outf = open(out,'w')
    c=0
    for line in inf:
	ll =line.strip().split("\t")
	if len(ll)<=5:
	    continue
	ID = ll[0]
	species = ll[1]
	factor = ll[2]
	folder = ll[3]
	QC = ll[4:]
		#folder = ll[2]
	xlsfile = folder+"/"+ID+"_peaks.xls"
	outfile = directory+ID+"_summits_cutoff_"+str(cutoff)+".bed"
	xlsselect(xlsfile,cutoff,outfile)
#print peakfile
#	cmd = "cp %s %s"%(summitfile,directory)
		#print cmd
#	c+=1
#	if c==40:
#	    inf.close()
#	    exit()

def xlsselect(xlsfile,cutoff,outfile):
    inf = open(xlsfile)
    outf = open(outfile,'w')
    n=1
    for line in inf:
	name = "MACS_peak_%s"%(str(n))
	n+=1
	if line.startswith("#") or line.strip()=="":
	    continue	
	ll = line.split()
	if ll[0]=="chr" :
	    continue
	if len(ll)==9:
	    peakname=name
	else:
	    peakname = ll[9]
	FC = ll[7]
	if float(FC)>=cutoff:
	    try:
	        newline = "\t".join([ll[0],str(int(ll[4])-150),str(int(ll[4])+150),peakname,ll[8]])+"\n"
	    except:
		print line
		exit()
	    outf.write(newline)

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog <-p deletepath> "
    description = """Step 1 for Chongzhi work, fetch all dnaseI peak from server3, using DCzoho sheet input."""

    optparser = OptionParser(version="%prog 0.72",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========core options=============
    optparser.add_option("-d","--directory",dest="directory",type="str",
                         help="")
    optparser.add_option("-t","--DNtable",dest="table",type="str",
                             help="")
#    optparser.add_option("-o","--out",dest="out",type="str",
#                                 help="")
    optparser.add_option("-c","--cutoff",dest="cutoff",type="int",
                                     help="")

    (options,args) = optparser.parse_args()

    directory = options.directory
    table = options.table
#    out = options.out
    cutoff = options.cutoff
    if not table:
        optparser.print_help()

        sys.exit(1)
    fetchpeak(table,directory,cutoff)



if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

