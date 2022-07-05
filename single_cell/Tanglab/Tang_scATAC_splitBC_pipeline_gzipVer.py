#!/usr/bin/env python2.7
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
import copy,time
#import numpy
#import twobitreader
import gzip
from Levenshtein import distance

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def single_barcode_adj(b,BClist):
    if b in BClist:
        return [b,"match"]
    else:
        dis_to_refer = []
        for referBC in BClist:
            dis_to_refer.append( distance(b,referBC) )
        sorted_dis_to_refer = sorted(dis_to_refer)
        if sorted_dis_to_refer[0] <= 1 and \
           sorted_dis_to_refer[1] - sorted_dis_to_refer[0] >= 1:# \
           #and dis_to_refer.index(sorted_dis_to_refer[0]) <= 1:
            return [BClist[dis_to_refer.index(sorted_dis_to_refer[0])],"adj"]
        else:
            return ["NA","rm"]
#    if b in P7_refer[:6]:
#        return b
#    else:
#        for referBC in P7_refer:
#            dis_to_refer.append( distance(b,referBC) )
#        sorted_dis_to_refer = sorted(dis_to_refer)
#        if sorted_dis_to_refer[0] <= 1 and \
#           sorted_dis_to_refer[1] - sorted_dis_to_refer[0] >= 1 :# \
#           #and dis_to_refer.index(sorted_dis_to_refer[0]) <= 3:
#            return P7_refer[dis_to_refer.index(sorted_dis_to_refer[0])]
#        else:
#            return "NA"   

def splitATACfq(inputname,metafile,folder,P5,P7,startpos):

    P5list = P5.split(",")
    P7list = P7.split(",")

    if len(set(map(len,P5list))) == 1:
        barcodeLen5 = len(P5list[0])
    else:
        print "P5 barcodes has different length"
        print "pipeline exit"
        sys.exit()

    if len(set(map(len,P7list))) == 1:
        barcodeLen7 = len(P7list[0])
    else:
        print "P7 barcodes has different length"
        print "pipeline exit"
        sys.exit()    

    ## readin meta data file, scan individual cell name corresponded to the fq file name
    inf = open(metafile)
    bc2cell = {}
    bc2cellEX = {}
    bc2cellname = {}
    for line in inf:
        ll = line.strip().split()
        if not ll[0] == inputname:
            #print ll
            continue
        BC = ll[3]
        p5 = BC.split("_")[0]
        p7 = BC.split("_")[1]
        if not p5 in P5list:
            print "P5 %s barcode is not in P5list"%(p5)
            print ll
            print "pipeline exit"
            sys.exit()
        if not p7 in P7list:
            print "P7 %s barcode is not in P7list"%(p7)
            print ll
            print "pipeline exit"
            sys.exit()

        if bc2cell.has_key(BC):
            print "barcode %s is not unique in fq: %s"%(BC, inputname)
            print "pipeline exit"
            sys.exit()
        else:
            bc2cell[BC] = [gzip.open(ll[5]+"_R1.fastq.gz",'wb'),gzip.open(ll[5]+"_R2.fastq.gz",'wb')]
            bc2cellEX[BC] = [gzip.open(ll[5]+"ex_R1.fastq.gz",'wb'),gzip.open(ll[5]+"ex_R2.fastq.gz",'wb')]
            bc2cellname[BC] = ll[5]
    inf.close()

    unmatchcell = [gzip.open(inputname+"_unmatch_R1.fastq.gz",'wb'),gzip.open(inputname+"_unmatch_R2.fastq.gz",'wb')]

    ## scan fastq, fetch barcode count list
    if not folder.endswith("/"):
        folder += "/"

    fq1 = "%s%s_R1.fastq.gz"%(folder,inputname)
    fq2 = "%s%s_R2.fastq.gz"%(folder,inputname)
    inf1 = gzip.open(fq1,"rb")
    inf2 = gzip.open(fq2,"rb")

    line_count = 0
    obs_barcode_reads_count = {}
    #print bc2cell.keys()
    while 1:
        line_count += 1
        p1 = inf1.readline()
        p2 = inf2.readline()

        if line_count%4 == 1:
            p1l1 = p1
            p2l1 = p2
        elif line_count%4 == 2:
            p1l2 = p1[startpos:]
            p2l2 = p2[startpos:]
            obs_p5 = p1[:barcodeLen5]
            obs_p7 = p2[:barcodeLen7] 
            #print "p1",p1,line_p5
            #print "p2",p2,line_p7               
        elif line_count%4 == 3:
            p1l3 = p1
            p2l3 = p2
        else:
            p1l4 = p1[startpos:]
            p2l4 = p2[startpos:]


            obs_barcode = obs_p5 + "_" + obs_p7

            ### add to BC summary
            if obs_barcode_reads_count.has_key(obs_barcode):
                use_barcode = obs_barcode_reads_count[obs_barcode][0]
            else:
                use_p5_item = single_barcode_adj(obs_p5,P5list)
                use_p7_item = single_barcode_adj(obs_p7,P7list)
                use_barcode = use_p5_item[0]+ "_" + use_p7_item[0]

                if bc2cellname.has_key(use_barcode):
                    individual_cellname = bc2cellname[use_barcode]
                else:
                    individual_cellname = "NA"
                if obs_barcode == use_barcode:
                    obs_barcode_reads_count[obs_barcode] = [use_barcode,"match",individual_cellname,0]
                elif not "NA" in use_barcode:
                    obs_barcode_reads_count[obs_barcode] = [use_barcode,"adj",individual_cellname+"ex",0]
                else:
                    obs_barcode_reads_count[obs_barcode] = [use_barcode,"rm",individual_cellname,0]

            obs_barcode_reads_count[obs_barcode][3] += 1 


            ### match individual cell name

            if bc2cell.has_key(use_barcode):
                if obs_barcode == use_barcode:
                    bc2cell[use_barcode][0].write(p1l1.split()[0]+":"+use_barcode+":"+obs_barcode+" "+p1l1.split()[1]+ "\n")
                    bc2cell[use_barcode][0].write(p1l2)
                    bc2cell[use_barcode][0].write(p1l3)
                    bc2cell[use_barcode][0].write(p1l4)
                    bc2cell[use_barcode][1].write(p2l1.split()[0]+":"+use_barcode+":"+obs_barcode+" "+p2l1.split()[1]+ "\n")
                    bc2cell[use_barcode][1].write(p2l2)
                    bc2cell[use_barcode][1].write(p2l3)
                    bc2cell[use_barcode][1].write(p2l4)
                    #print 1
                else:
                    bc2cellEX[use_barcode][0].write(p1l1.split()[0]+":"+use_barcode+":"+obs_barcode+" "+p1l1.split()[1]+ "\n")
                    bc2cellEX[use_barcode][0].write(p1l2)
                    bc2cellEX[use_barcode][0].write(p1l3)
                    bc2cellEX[use_barcode][0].write(p1l4)
                    bc2cellEX[use_barcode][1].write(p2l1.split()[0]+":"+use_barcode+":"+obs_barcode+" "+p2l1.split()[1]+ "\n")
                    bc2cellEX[use_barcode][1].write(p2l2)
                    bc2cellEX[use_barcode][1].write(p2l3)
                    bc2cellEX[use_barcode][1].write(p2l4)
            else:
                unmatchcell[0].write(p1l1.split()[0]+":"+use_barcode+":"+obs_barcode+" "+p1l1.split()[1]+ "\n")
                unmatchcell[0].write(p1l2)
                unmatchcell[0].write(p1l3)
                unmatchcell[0].write(p1l4)
                unmatchcell[1].write(p2l1.split()[0]+":"+use_barcode+":"+obs_barcode+" "+p2l1.split()[1]+ "\n")
                unmatchcell[1].write(p2l2)
                unmatchcell[1].write(p2l3)
                unmatchcell[1].write(p2l4)

        #if line_count == 1000000:
        #    break
        if p1.strip() == "" and p2.strip() == "":
            break
    inf1.close()
    inf2.close()
    unmatchcell[0].close()
    unmatchcell[1].close()

    for BC in bc2cell.keys():
        #print BC,len(bc2cell[BC]),len(bc2cellEX[BC])
        bc2cell[BC][0].close()
        bc2cell[BC][1].close()
        bc2cellEX[BC][0].close()
        bc2cellEX[BC][1].close()

    #print obs_barcode_reads_count
    #for i in obs_barcode_reads_count.iteritems():
    #    print i
    outf = open("%s_BCcount.txt"%inputname,'w')

    for obsBC,features in sorted(obs_barcode_reads_count.iteritems(),key=lambda x:x[1][3],reverse=True):
        outf.write("\t".join(map(str,[obsBC]+features))+"\n")
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


    optparser.add_option("-m",dest="metafile",type="str",
                         help="file for meta data")
    optparser.add_option("-i",dest="inputname",type="str",
                         help="PE fastq name is: inputname_R1.fastq, inputname_R2.fastq. name should also in metadata file, col1")
    optparser.add_option("--p5",dest="p5",type="str",default="TATAGCCT,ATAGAGGC,CCTATCCT,GGCTCTGA,AGGCGAAG,TAATCTTA,CAGGACGT,GTACTGAC",
                         help="list of P5barcodes")
    optparser.add_option("--p7",dest="p7",type="str",default="CGAGTAAT,TCTCCGGA,AATGAGCG,GGAATCTC,TTCTGAAT,ACGAATTC,AGCTTCAG,GCGCATTA,CATAGCCG,TTCGCGGA,GCGCGAGA,CTATCGCT",
                         help="list of P7barcodes")
    optparser.add_option("--startpos",dest="startpos",type="int",default=27,
                         help="")
    optparser.add_option("-f","--folder",dest="folder",type="str",default="/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/fastq/06_20181122_mESC_CAT4/combined_ATAC/CAT4_P1",
                         help="folder for fastq files")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputname = options.inputname
    if not inputname:
        optparser.print_help()
        sys.exit(1)
    splitATACfq(inputname,options.metafile,options.folder,options.p5,options.p7,options.startpos)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



