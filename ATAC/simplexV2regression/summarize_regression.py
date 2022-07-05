import os,sys

def fetch_rf_cv(name):
    inf = open(name)
    line = inf.readline()
    line = inf.readline()
    #for line in inf:
    #    if line.startswith("rf"):
    cv = line.split()[-2]#CVdict[name][3] = line.split()[2]
    inf.close()
    return cv
CVrf = {}

folder = 'ATAC_simplexV2mat_ATACbias/'

for f in sorted(os.listdir(folder)):
    if f.endswith("_simplexV2output.bed"):
        name = f.split("_simplexV2output.bed")[0]
        CVrf[name] = ["NA"]*10
        ## ATAC_ATAC,ATAC_fakeDNase,ATAC_fakeRD,
        ## ATAC_permuteSeq,ATAC_randomSeq
        ## DNase_DNase,DNase_fakeATAC
print CVrf.keys()
folders = [ 'ATAC_simplexV2mat_ATACbias/',
            'ATAC_simplexV2mat_fakeDNasebias/',
            'ATAC_simplexV2mat_fakeRDbias/',
            'ATAC_simplexV2mat_permuteSequence/',
            'ATAC_simplexV2mat_randomSequence/',
            'DNase_simplexV2mat_DNasebias/',
            'DNase_simplexV2mat_fakeATACbias/',
            '../simplexVer2_regression/']

for i in range(len(folders)):
    folder = folders[i]
    for f in sorted(os.listdir(folder)):
        if f.endswith("plus_rf_CVpower.txt"):
            name = f.split("_plus_rf_CVpower.txt")[0]
            if i == 7:
                if not name.endswith("_10000"):
                    continue

            if os.path.isfile(folder + name+"_plus_rf_CVpower.txt"):
                this_cv = fetch_rf_cv(folder + name+"_plus_rf_CVpower.txt")
                if not name.startswith("GM12878"):
                    name  = "GM12878_" + name
                CVrf[name.split("_10000")[0]][i] = this_cv


folder = 'totalreads_regression/'
for f in sorted(os.listdir(folder)):
    if f.endswith("ATACsig_rf_CVpower.txt"):
        name = f.split("_ATACsig_rf_CVpower.txt")[0]
        if os.path.isfile(folder + name+"_ATACsig_rf_CVpower.txt"):
            this_cv = fetch_rf_cv(folder + name+"_ATACsig_rf_CVpower.txt")
            if not name.startswith("GM12878"):
                name  = "GM12878_" + name
            CVrf[name.split("_10000")[0]][8] = this_cv
    elif f.endswith("DNasesig_rf_CVpower.txt"):
        name = f.split("_DNasesig_rf_CVpower.txt")[0]
        if os.path.isfile(folder + name+"_DNasesig_rf_CVpower.txt"):
            this_cv = fetch_rf_cv(folder + name+"_DNasesig_rf_CVpower.txt")
            if not name.startswith("GM12878"):
                name  = "GM12878_" + name
            CVrf[name.split("_10000")[0]][9] = this_cv



outf = open("motif10k_rfCVpower_summary.txt",'w')
newll = ["name","ATAC_ATAC","ATAC_fakeDNase","ATAC_fakeRD",
            "ATAC_permuteSeq","ATAC_randomSeq",
            "DNase_DNase","DNase_fakeATAC","ATAColdnoGC",
            "ATACsig","DNasesig"]
outf.write("\t".join(newll)+"\n")
for motif in sorted(CVrf.keys()):
    newll = [motif] + CVrf[motif]
    outf.write("\t".join(newll)+"\n")
outf.close()



#for f in sorted(os.listdir(".")):
#    if f.endswith("plus_rf_CVpower.txt"):
#        name = f.split("_plus_rf_CVpower.txt")[0]
#        if os.path.isfile(name+"_plus_rf_CVpower.txt") \
#            and os.path.isfile(name+"_plus_lm_CVpower.txt"):
#            CVdict[name] = [0,0,0,0]
#            inf = open(name+"_plus_lm_CVpower.txt")
#            for line in inf:
#                if line.startswith("rawlm"):
#                    CVdict[name][0] = line.split()[2]
#                elif line.startswith("step"):
#                    CVdict[name][1] = line.split()[2]
#                elif line.startswith("lasso"):
#                    CVdict[name][2] = line.split()[2]
#            inf.close()
#            inf = open(name+"_plus_rf_CVpower.txt")
#            for line in inf:
#                if line.startswith("rf"):
#                    CVdict[name][3] = line.split()[2]
#            inf.close()
#outf = open("motif10k_CVpower_summary.txt",'w')
#newll = ["name","rawlm","stepwise","lasso","randomforest"]
#outf.write("\t".join(newll)+"\n")
#for motif in sorted(CVdict.keys()):
#    newll = [motif] + CVdict[motif]
#    outf.write("\t".join(newll)+"\n")
#outf.close()









