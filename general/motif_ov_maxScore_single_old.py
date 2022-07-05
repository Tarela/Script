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
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit() 
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

def add_ov_motif_score(inputpeak,motif_folder,motif_names,outname):
    ### overlap
    inputdict = {}
    inf = open(inputpeak)
    for line in inf:
        ll = line.split()
        inputdict[ll[3]] = ll
    inf.close()
    if not motif_folder.endswith('/'):
        motif_folder += '/'
    
    allmotif = motif_names.split(",")
    for motif in allmotif:
        cmd = 'intersectBed -a %s -b %s -wao > %s'%(inputpeak,motif_folder + motif + '.bed',outname+'.tmp')
        sp(cmd)
        inf = open(outname+'.tmp')
        this_dict = {}
        for line in inf:
            ll = line.split()
            if not this_dict.has_key(ll[3]):
                this_dict[ll[3]] = []
            if ll[-3] == '-1':
                ll[-3] = '0'
            this_dict[ll[3]].append(float(ll[-3]))
        inf.close()
        for peak in inputdict.keys():
            if not this_dict.has_key(peak):
                print motif,peak
            inputdict[peak].append(max(this_dict[peak]))
            
    outf = open(outname,'w')
    newll = ['chrm','start','end','name','disLast'] + allmotif
    outf.write("\t".join(newll)+"\n")
    for peak in sorted(inputdict.keys()):
        outf.write("\t".join(map(str,inputdict[peak]))+"\n")
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
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="")
    optparser.add_option("-f",'--motiffolder',dest="motif_folder",type="str",default = "/mnt/Storage/home/huse/Data/Motif/site/jaspar_v_zv9_uniq/",
                         help="motif folder, default: /mnt/Storage/home/huse/Data/Motif/site/jaspar_v_zv9_uniq/")
    optparser.add_option("-m",'--motifnames',dest="motif_names",type="str",default = "",
                         help="motif name, separate by comma, eg: MC001_AR,MC002_CTCF")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile

    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    add_ov_motif_score(inputfile,options.motif_folder,options.motif_names,options.output)
    #get_signal(inputfile,output,signalbw)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


