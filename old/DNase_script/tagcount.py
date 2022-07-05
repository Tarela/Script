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
import subprocess

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
def fetchpeak(table,directory,cmdfile,binfile):
    if not directory.endswith("/"):
    	directory+="/"
#    countfile
    inf = open(table,'r')
    outf = open(cmdfile,'w')
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
	allfile = os.listdir(folder)
	bamfiles = []	
	for f in allfile:
	    if f.startswith(str(ID)+"_rep") and f.endswith("_treat.bam"):
		bamfiles.append(folder+"/"+f)
#		print f
#	print bamfiles
	if len(bamfiles)==1:
	    cmd = "cp %s %s"%(bamfiles[0],directory+ID+".bam")
#	    print cmd
	elif len(bamfiles)>1:
	    cmd = "samtools merge %s %s"%(directory+ID+".bam"," ".join(bamfiles))
	else:
	    print ID+" bug in bamfile"
	    continue
#        os.system(cmd)
	cmd2 = "coverageBed -abam %s -b %s > %s"%(directory+ID+".bam",binfile,directory+ID+"_cover.bed") 
#	os.system(cmd2)
	cmd3 = "bamToBed -i %s > %s"%(directory+ID+".bam",directory+ID+".bed")
	cmd4 = "rm %s"%(directory+ID+".bam")
	cmd5 = "wc -l %s > %s "%(directory+ID+".bed",directory+ID+"_totalreads.txt")
	cmd6 = "rm %s"%(directory+ID+".bed")
#	os.system(cmd3)
	outf.write(cmd+"\n"+cmd2+"\n"+cmd3+"\n"+cmd4+"\n"+cmd5+"\n"+cmd6+"\n")
#	cmd4 = "wc -l %s"%(directory+ID+".bed")
#	totalreads = sp(cmd4)[0].strip().split()[0]
#	outf.write(ID+"\t"+totalreads+"\n")
    inf.close()
    outf.close()
#	print allfile
#	bamfile = folder+"/"+ID+"_peaks.xls"
#	outfile = directory+ID+"_summits_cutoff_"+str(cutoff)+".bed"
#	xlsselect(xlsfile,cutoff,outfile)
#print peakfile
#	cmd = "cp %s %s"%(summitfile,directory)
		#print cmd
#	c+=1
#	if c==10:
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
    optparser.add_option("-o","--out",dest="out",type="str",
                                 help="")
    optparser.add_option("-b","--binfile",dest="binfile",type="str",
                                     help="")

    (options,args) = optparser.parse_args()

    directory = options.directory
    table = options.table
    out = options.out
    binfile = options.binfile
    if not table:
        optparser.print_help()

        sys.exit(1)
    fetchpeak(table,directory,out,binfile)



if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

