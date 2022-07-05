#!/usr/bin/env python
'''
Created on XXXX-XX-XX

@author: Tarela
'''
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging,time
import string
import re 
import urllib
from ftplib import FTP
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def addSRA(infile,outname,nGSMcol):

    inf = open(infile)
    outmeta = open(outname+"_meta.txt",'w')
    outcmd = open(outname+'_cmd.sh','w')
    
    ftp_connection = FTP('ftp-trace.ncbi.nlm.nih.gov')
    ftp_connection.login()
    F = ftp_connection
    
    for line in inf:
        ll = line.split()
        gsm = ll[max(0,int(nGSMcol)-1)].strip()
        if not gsm.startswith("GSM"):
            print "GSMID is not startswith GSM",ll
        weblink = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s'%(gsm)
        webpage = urllib.urlopen(weblink).read().split("\n")
        SRXID = ""
        for line in webpage:
            if 'www.ncbi.nlm.nih.gov/sra?term=SRX' in line:
                SRXstart = line.find('>SRX')
                SRXID = line[(SRXstart+1):].split("<")[0]
                #SRXstart = line.find('www.ncbi.nlm.nih.gov/sra?term=SRX')
                #SRXID = line[(SRXstart+37):(SRXstart+46)]
                #for srxlength in range(8,11):
                #    if line[(SRXstart+srxlength)] in ['"','<']:
                #        SRXID = line[SRXstart:(SRXstart+srxlength)]
                #        break
                break
        if SRXID == "":
            print gsm,"noSRXID"
            continue
        elif not SRXID.startswith('SRX'):
            print gsm,SRXID,'errorID'
            continue    
        else:
            try:
                F.cwd('/sra/sra-instant/reads/ByExp/sra/SRX/%s/%s'%(SRXID[0:6],SRXID))
            except:
                print SRXID,'noSRRID'
                continue
        SRRIDs = F.nlst()
        for srrid in SRRIDs:
            outcmd.write(" ".join(["./downmap",SRXID,srrid])+"\n")
        newll = ll + [SRXID,",".join(SRRIDs)]
        outmeta.write("\t".join(newll)+"\n")
    outmeta.close()
    outcmd.close()
        	


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """add sequence feature following each line,(ext<-center->ext), newll  = ll+[localMaxCpG,aveCpG,GCcontent]"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",
                         help="")
    optparser.add_option("-o","--outname",dest="outname",type="str",default="GSM2SRR",
                         help="name of output files, you get two output files: outname_meta.txt, outname_cmd.sh, default is GSM2SRR")
    optparser.add_option("-n","--nGSMcol",dest="nGSMcol",type="int",default=1,
                         help="column number of GSMID, default is 1 (first column)")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    outname = options.outname
    nGSMcol = options.nGSMcol
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    addSRA(inputfile,outname,nGSMcol)
    #get_signal(inputfile,output,signalbw,extdown,extup,number)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
