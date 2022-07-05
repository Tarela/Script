import subprocess
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
import os
inf = open("motif10k_rfCVpower_summary.txt")
for line in inf:
    ll = line.split()
    name = ll[0]
    rfcv = ll[1]
    if name == "name":
        continue
    if ll[1] == "NA":
        continue
    print "\nATAC\t"+name+"\t"+rfcv
    #cmd = "awk '{if ($5 > 10000) print $0}' ../oneDim_GM12878/GM12878_%s_ChIPpeak_ATACpeak_mappable.bed | wc -l "%(name)
    #motifnum = sp(cmd)[0].split()[0]
    #print "\n"+name+"\t"+rfcv+"\t"+motifnum
    #os.system(cmd)
    cmd = "sort -k 2,2gr ATAC_simplexV2mat_ATACbias/%s_plus_rf_importance.txt |head"%(name)
    os.system(cmd)
    print "\nDNase\t"+name+"\t"+rfcv
    #cmd = "awk '{if ($5 > 10000) print $0}' ../oneDim_GM12878/GM12878_%s_ChIPpeak_ATACpeak_mappable.bed | wc -l "%(name)
    #motifnum = sp(cmd)[0].split()[0]
    #print "\n"+name+"\t"+rfcv+"\t"+motifnum
    #os.system(cmd)
    cmd = "sort -k 2,2gr DNase_simplexV2mat_DNasebias/%s_plus_rf_importance.txt |head"%(name)
    os.system(cmd)

