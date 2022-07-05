#!/usr/bin/env python
#Time-stamp:<Tarela>
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

def readTag(tagfile,t):
    tagdict = {}
    inf = open(tagfile)
    for line in inf:
        ll = line.split()
        if len(ll) == 0:
            continue
        if not tagdict.has_key(ll[0]):
            tagdict[ll[0]] = []
        if len(ll)>=6 and ll[5] == '-':
            if t == 'inte':
                tagdict[ll[0]].append(int(ll[2])-1)
            else :
                tagdict[ll[0]].append(float(ll[2])-1)
        else:
            if t == 'inte':
                tagdict[ll[0]].append(int(ll[1]))
            else:
                tagdict[ll[0]].append(float(ll[1]))
    inf.close()    
    return tagdict

def getdistance(taga,tagb,out,maxdis):
    maxdis = abs(maxdis)
    ta = readTag(taga,'inte')
    tb = readTag(tagb,'float')
    c=0
    alldis = {}
    print 'start'
    for chrm in ta.keys():
        if not tb.has_key(chrm):
            continue
        v = sorted(ta[chrm] + tb[chrm])
     
        for i in range(len(v)):
            c+=1
            if c%100000 ==0:
                print c
            if type(v[i]) == int:
                j = i
                while j>0:
                    j -= 1
                    if type(v[j]) != float:
                        continue
                    dis = v[i]-v[j]
                    if dis <= maxdis:
                        if not alldis.has_key(dis):
                            alldis[dis] = 0
                        alldis[dis] += 1
                    else:
                        break
                j=i
                while j < len(v)-1:
                    j += 1                
                    if type(v[j]) != float:
                        continue
                    dis = v[i]-v[j]
                    
                    if dis >= -maxdis:
                        if not alldis.has_key(dis):
                            alldis[dis] = 0
                        alldis[dis] += 1  
                    else:
                        break
    outf = open(out,'w')
    for i in sorted(alldis.keys()):
        newll =  [i,alldis[i]]
        outf.write("\t".join(map(str,newll))+"\n")
    outf.close()
                    
        #for i in range(len(va)):
    

# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-a","--taga",dest="taga",type="str",
                         help="")
    optparser.add_option("-b","--tagb",dest="tagb",type="str",
                         help="")
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")

#========minor options=============
    optparser.add_option("-m","--maxdis",dest="maxdis",type="int",default=300,
                         help="maxdistance from tag1 to tag2, dafault is 300 , means calculate tag2 locate in +-300bp from tag1")


    (options,args) = optparser.parse_args()

    taga = options.taga
    tagb = options.tagb
    out = options.out
    maxdis = options.maxdis
    if not taga:
        optparser.print_help()
        sys.exit(1)
    
    getdistance(taga,tagb,out,maxdis)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

