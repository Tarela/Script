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
import random

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
import subprocess
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
def sample(inputfile,outputfile,percentage,number):
    if percentage == -1 :
        total = int(sp('wc -l %s'%(inputfile))[0].split()[0])
        percentage = number*100.0/total
    inf = open(inputfile)
    outf = open(outputfile,'w')
    for line in inf:
        a = random.randint(0,10000)
        if a<= percentage*100:
            outf.write(line)
    inf.close()
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
    optparser.add_option("-p","--percentage",dest="percentage",type="int",default = -1,
                         help="% of reads, 10 = 10%, can't be set with -n")
    optparser.add_option("-n","--number",dest="number",type="int",default = -1,
                         help="number of reads , can't be set with -p ")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    percentage = options.percentage
    number = options.number
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    if percentage ==  -1 :
        if number == -1 :
            print 'none of -p and -n are set, only one of them could exist'
            optparser.print_help()
            sys.exit(1)
            
    else : 
        if number != -1:
            print 'both -p and -n are set ,only one of them could exist '
            optparser.print_help()
            sys.exit(1)

    sample(inputfile,output,percentage,number)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

