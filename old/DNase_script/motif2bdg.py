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
import math
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def convert(inputbed, outputbdg,uselog):
    input_dict = {}
    inf = open(inputbed)
    for line in inf:
        ll = line.split()
        if input_dict.has_key(ll[0]):
            if uselog ==1 :
                input_dict[ll[0]].append([int(ll[1]),int(ll[2]),math.log10(float(ll[4]))])
            else:
                input_dict[ll[0]].append([int(ll[1]),int(ll[2]),(float(ll[4]))])
        else:
            if uselog==1:
                input_dict[ll[0]]=[[int(ll[1]),int(ll[2]),math.log10(float(ll[4]))]]
            else:
                input_dict[ll[0]]=[[int(ll[1]),int(ll[2]),(float(ll[4]))]]
    inf.close()
    outf = open(outputbdg,'w')
    for chrm in input_dict.keys():
        l = sorted(input_dict[chrm],key=lambda x:x[0])
        for i in range(1,len(l)):
            if l[i][0]>=l[i-1][1]:
                ll = [chrm]+l[i-1]
                outf.write("\t".join(map(str,ll))+"\n")
            else:
                ll = [chrm,l[i-1][0],l[i][0],l[i-1][2]]
                outf.write("\t".join(map(str,ll))+"\n")
                #ll = [chrm,l[i][0],l[i-1][1],l[i-1][2]+l[i][2]]
                #outf.write("\t".join(map(str,ll))+"\n")
                #l[i][0]=l[i-1][1]
    outf.close()

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog <-p deletepath> "
    description = """delete all unuse pipeline result in given path"""

    optparser = OptionParser(version="%prog 0.72",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========core options=============
    optparser.add_option("-i","--inputbed",dest="inputbed",type="str",
                         help="")
    optparser.add_option("-o","--outputbdg",dest="outputbdg",type="str",
                         help="")
    optparser.add_option("--log",dest="uselog",type="int",default=1,
                         help="1 log , 0 linear, default 1")

    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    outputbdg = options.outputbdg
    uselog = options.uselog
    if not inputbed or not outputbdg:
        optparser.print_help()

        sys.exit(1)
    convert(inputbed, outputbdg,uselog)



if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

