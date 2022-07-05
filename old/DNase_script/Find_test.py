'''
Created on 2013-01-02

@author: Shawn
'''
#!/usr/bin/env python
#Time-stamp:<2013-01-02 Shawn>
"""
Description: >.<

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging
import string
import random
import time
import scipy
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def read_cleavage(inputfile,nostrand):
    cleavage = {}
    location = {}
    inf = open(inputfile)
    i=0
    for line in inf:
        i += 1
        name = "Chipsite_"+str(i)
        ll = line.split()
        ll[3]=name
        location[name]=ll
        if not nostrand:
            plus = map(float,ll[6].split(","))
            minus = map(float,ll[7].split(","))
            total = list(map(lambda x: x[0]+x[1], zip(plus, minus))) 
            cleavage[name]=total
        else:
            cleavage[name]=map(float,ll[6].split(","))
    return cleavage,location

def caculate_footprint(digital_list,central_max,central_min,flanking_max,flanking_min,site,cutoff):
    footprint = []
    last_start = "NA"
    last_end = "NA"
    last_value = "NA"
    for i in range(len(digital_list)): ### i is current position -> right bound of central region
        if i < flanking_min or len(digital_list) - flanking_min - central_min < i :
            continue
        for central_length in range(central_min,central_max):
            for flanking_length in range(flanking_min,flanking_max):
                l_end = i
                c_start = i
                if len(digital_list) - flanking_min > i + central_length:
                    c_end = i + central_length
                else:
                    c_end = len(digital_list) - flanking_min
                r_start = c_end
                r_end = r_start + flanking_length
                if r_end > len(digital_list):
                    r_end = len(digital_list)
                    l_start = i - ( r_end - r_start )
                elif i <= flanking_length:
                    l_start = 0
                    r_end = r_start + i
                else:
                    l_start = i - flanking_length

                L = scipy.mean(digital_list[l_start:l_end])
                C = scipy.mean(digital_list[c_start:c_end])
                R = scipy.mean(digital_list[r_start:r_end])
                if L <= C or R <= C:
                    continue
                else:
                    FOS = (C+1)/L + (C+1)/R
                if last_start == "NA" :
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
                elif c_start >= last_end:
                    if last_value < cutoff :
                        footprint.append([last_start,last_end,last_value,site])
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
                elif FOS < last_value:
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
    if last_value < cutoff:
        footprint.append([last_start,last_end,last_value,site])
    return footprint
def caculate_footprint_value(digital_list,central_max,central_min,flanking_max,flanking_min):
    footprint = []
    last_start = "NA"
    last_end = "NA"
    last_value = "NA"
    for i in range(len(digital_list)): ### i is current position -> right bound of central region
        if i < flanking_min or len(digital_list) - flanking_min - central_min < i :
            continue
        for central_length in range(central_min,central_max):
            for flanking_length in range(flanking_min,flanking_max):
                l_end = i
                c_start = i
                if len(digital_list) - flanking_min > i + central_length:
                    c_end = i + central_length
                else:
                    c_end = len(digital_list) - flanking_min
                r_start = c_end
                r_end = r_start + flanking_length
                if r_end > len(digital_list):
                    r_end = len(digital_list)
                    l_start = i - ( r_end - r_start )
                elif i <= flanking_length:
                    l_start = 0
                    r_end = r_start + i
                else:
                    l_start = i - flanking_length

                L = scipy.mean(digital_list[l_start:l_end])
                C = scipy.mean(digital_list[c_start:c_end])
                R = scipy.mean(digital_list[r_start:r_end])
                if L <= C or R <= C:
                    continue
                else:
                    FOS = (C+1)/L + (C+1)/R
                if last_start == "NA" :
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
                elif c_start >= last_end:
                    footprint.append(last_value)
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
                elif FOS < last_value:
                    last_start = c_start
                    last_end = c_end
                    last_value = FOS
    footprint.append(last_value)
    return footprint


def summary(cleavage,location,central_max,central_min,flanking_max,flanking_min,cutoff):
    footprints = []    #footprints =[[start,end,value,sitename]]
    shuffle_value = []
    for site in cleavage:
        t=time.time()
        digital = cleavage[site]
        footprints.extend(caculate_footprint(digital,central_max,central_min,flanking_max,flanking_min,site,cutoff))
#        print time.time()-t
    #    t=time.time()
    #    random.shuffle(digital)
    #    shuffle_value.extend(caculate_footprint_value(digital,central_max,central_min,flanking_max,flanking_min))
#        print time.time()-t
    footprints.sort(key = lambda x:x[2])
    #print footprints
    #shuffle_value.sort()
    #print shuffle_value[1:5]
    #i=0
    #j=0
    #while i < len(footprints):
     #   print i,j,footprints[i]
    #    if footprints[i][2] < shuffle_value[j]:
    #        if j==0:
    #            footprints[i].append(0.0)
    #            footprints[i].append(i)
    #            footprints[i].append(j)
    #        else:
    #            footprints[i].append((j+0.0001)/(i+0.0001))
    #            footprints[i].append(i)
    #            footprints[i].append(j)
    #        i += 1
    #    else:
    #        if j < len(shuffle_value)-1:
    #            j += 1
    #        else:
    #            if j==0:
    #                footprints[i].append(0.0)
    #                footprints[i].append(i)
    #                footprints[i].append(j)
    #            else:
    #                footprints[i].append((j+0.0001)/(i+0.0001))
    #                footprints[i].append(i)
    #                footprints[i].append(j)
    #            i += 1
    #print sorted(footprints,key=lambda x:x[2])
    outf = open('a1k.r','w')
    l=[]
    for i in footprints:
        l.append(i[2])
    outf.write("a<-c("+str(l)[1:-1]+")\n")
    outf.write("plot(dentisy(a))\ndev.off()\n")
    outf.close()
        
    #for site   


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog <-i standard input> <-o output> [options] "
    description = """Description : find footprint in digital dnase cleavage data, the input data should have more than 8 column , 1~5 is macs peak output ,6th column is the mean dnase cleavage density, 7,8th is digital dnase cleavage (+/-) , if dnase cleavge do not sep strand, specify the --nostrand param to let the script read only column 6 for total cleavage"""

    optparser = OptionParser(version="%prog 0.72",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--inputfile",dest="inputfile",type="str",
                         help="standard input file")
    optparser.add_option("--nostrand",dest="nostrand",action='store_true',
                         help="if set, script only read column 7 for digital dnase cleavage(no strand) ; else script will sum up the 7th and 8th value for digital dnase cleavage(+/- strand)")
    optparser.add_option("-o","--output",dest="output",type="str",
                         help="output file name")

#========minor options============
    optparser.add_option("--c_max",dest="c_max",type="int",default = 25,
                         help="max length of central region")
    optparser.add_option("--c_min",dest="c_min",type="int",default = 6,
                         help="min length of central region")
    optparser.add_option("--f_max",dest="f_max",type="int",default = 10,
                         help="max length of flanking region")
    optparser.add_option("--f_min",dest="f_min",type="int",default = 3,
                         help="min length of flanking region")
    optparser.add_option("-c","--cutoff",dest="cutoff",type="int",default = 3,
                         help="output file name")
 
    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    nostrand = options.nostrand
    cutoff = options.cutoff
    c_max = options.c_max
    c_min = options.c_min
    f_max = options.f_max
    f_min = options.f_min
    if not inputfile or not output:
        optparser.print_help()
        sys.exit(1)

    cleavage,location = read_cleavage(inputfile,nostrand)
    summary(cleavage,location,c_max,c_min,f_max,f_min,cutoff)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

