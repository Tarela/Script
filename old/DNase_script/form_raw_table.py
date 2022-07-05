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

def form(table,directory,unionfile,rpkm,output):
    tagcount_dict={}
    readnum_dict = {}
    peakname = []
    peaklength = []
    inf = open(unionfile)
    for line in inf:
	ll = line.split()
	length = int(ll[2])-int(ll[1])
	name = "_".join(ll)
	peakname.append(name)
	peaklength.append(length)
    inf.close()
    if not directory.endswith("/"):
    	directory+="/"
    os.chdir(directory)
    inf = open(table,'r')
#    outf = open(out,'w')
    c=0
    for line in inf:
	ll =line.split()
	ID = ll[1].strip().split("/")[-1].split(".")[0]
	reads = ll[0]
	readnum_dict[ID]=reads
	tcfile = ID+"_cover.bed"
	tagcount_dict[ID]=get_tc(tcfile)
    inf.close()
    outf = open(output,'w')
    outf.write("\t".join(sorted(tagcount_dict.keys()))+"\n")
    for i in range(len(peakname)):
	pn = peakname[i]
	pl = peaklength[i]
	newll = [pn]
	for ids in sorted(tagcount_dict.keys()):
	    rn = int(readnum_dict[ids])
	    if rpkm:
		value = str(float(tagcount_dict[ids][i])*1000*100000/pl/rn)
	    else:
		value = str(tagcount_dict[ids][i])
	    newll.append(value)
	newline = "\t".join(newll)+"\n"
	outf.write(newline)
    outf.close()

def get_tc(tcfile):
    tcs=[]
    inf = open(tcfile)
    for line in inf:
	ll=line.split()
	tc=ll[-4]
	tcs.append(tc)
    inf.close()
    return tcs
#	species = ll[1]
#	factor = ll[2]
#	folder = ll[3]
#	QC = ll[4:]
#		#folder = ll[2]
#	xlsfile = folder+"/"+ID+"_peaks.xls"
#	outfile = directory+ID+"_summits_cutoff_"+str(cutoff)+".bed"
#	xlsselect(xlsfile,cutoff,outfile)
	
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
                         help="directory for all necessary files")
    optparser.add_option("-u","--union",dest="union",type="str",
                         help="union site file, given peak length and peakname ,name = chr+start+end")

    optparser.add_option("-t","--tagcounttable",dest="table",type="str",
                             help="a table show all sample and its tagcount, file should be in directory given by -d, only inpu file name, no path ~")
    optparser.add_option("--rpkm",dest="rpkm",action='store_true',
                         help="if --rpkm is input, we will caculate rpkm for each sample&site")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="output formal table ")

    (options,args) = optparser.parse_args()

    directory = options.directory
    table = options.table
    rpkm = options.rpkm
    union = options.union
    out = options.out
    if not table:
        optparser.print_help()

        sys.exit(1)
    form(table,directory,union,rpkm,out)



if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

