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
def sample(inputfile1,inputfile2,outputfile1,outputfile2,percentage,number):
#    random.seed(1228)
    if percentage == -1 :
        total = int(sp('wc -l %s'%(inputfile1))[0].split()[0])*1.0/4
        percentage = number*1000.0/total
        print percentage
    inf1 = open(inputfile1)
    inf2 = open(inputfile2)
    outf1 = open(outputfile1,'w')
    outf2 = open(outputfile2,'w')

    count= 0
    while 1:
        r1_p1=inf1.readline()
        if r1_p1.strip() != "":
            r2_p1=inf2.readline()
            r1_p2=inf1.readline()
            r2_p2=inf2.readline()
            r1_p3=inf1.readline()
            r2_p3=inf2.readline()
            r1_p4=inf1.readline()
            r2_p4=inf2.readline()
            count += 1
            a = random.randint(0,1000)
            if a<= percentage:
#                outf.write(line)
                outf1.write(r1_p1)
                outf1.write(r1_p2)
                outf1.write(r1_p3)
                outf1.write(r1_p4)
                outf2.write(r2_p1)
                outf2.write(r2_p2)
                outf2.write(r2_p3)
                outf2.write(r2_p4)

        else:
            break

    inf1.close()
    inf2.close()
    outf1.close()
    outf2.close()





# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============

    optparser.add_option("--inf1",dest="inputfile1",type="str",default = "",
                         help="")
    optparser.add_option("--inf2",dest="inputfile2",type="str",default = "",
                         help="")
    optparser.add_option("--outf1",dest="output1",type="str",default = "",
                         help="")
    optparser.add_option("--outf2",dest="output2",type="str",default = "",
                         help="")
    optparser.add_option("-p","--percentage",dest="percentage",type="float",default = -1,
                         help="% of reads, input 0.1 for 10%, can't be set with -n")
    optparser.add_option("-n","--number",dest="number",type="int",default = -1,
                         help="number of reads , can't be set with -p ")
#========minor options=============

    (options,args) = optparser.parse_args()

    percentage = options.percentage
    number = options.number
    if not options.inputfile1:
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

    sample(options.inputfile1,options.inputfile2,options.output1,options.output2,percentage,number)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

