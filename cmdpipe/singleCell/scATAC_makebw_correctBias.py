import os

import subprocess
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
#outf = open('line_count.txt','w')
folder = '/scratch/sh8tv/Project/scATAC/Data/scATAC/Nature2015_GSE65360/combineRep/'
for f in sorted(os.listdir(folder)):
    readcount = int(sp('wc -l %s'%(folder+f))[0].split()[0])
    if f.endswith('.bed') and os.path.isfile(folder + f) :
        name = f[:-4]
        cmd1 = """awk '{OFS="\\t";print $1,$2,$2+1,$4,".","+";}' %s > %s"""%(folder+f,name+"_plus.bed")
        cmd2 = """awk '{OFS="\\t";print $1,$3-1,$3,$4,".","-";}' %s > %s"""%(folder+f,name+"_minus.bed")
        cmd3 = """macs2 pileup -i %s -o %s --extsize 1 -f BED"""%(name+"_plus.bed", name+"_plus.bdg")
        cmd4 = """macs2 pileup -i %s -o %s --extsize 1 -f BED"""%(name+"_minus.bed", name+"_minus.bdg")
        cmd5 = """bdg2bw %s /scratch/sh8tv/Data/Genome/hg38/hg38_clean.len"""%(name+"_plus.bdg")
        cmd6 = """bdg2bw %s /scratch/sh8tv/Data/Genome/hg38/hg38_clean.len"""%(name+"_minus.bdg")
        cmd7 = """python   /scratch/sh8tv/Script/ATAC/correct_bias_bdg.py  -i %s -o %s -s + -b /scratch/sh8tv/Project/scATAC/Result/diffmer_var_comparison/diffmerMat/flanking/humanTcellNaked_ATACPE_enc4flank.txt -t fxr --genome /scratch/sh8tv/Data/Genome/hg38/hg38.2bit"""%(name+"_plus.bdg", name+"_plus")
        cmd8 = """python   /scratch/sh8tv/Script/ATAC/correct_bias_bdg.py  -i %s -o %s -s - -b /scratch/sh8tv/Project/scATAC/Result/diffmer_var_comparison/diffmerMat/flanking/humanTcellNaked_ATACPE_enc4flank.txt -t fxr --genome /scratch/sh8tv/Data/Genome/hg38/hg38.2bit"""%(name+"_minus.bdg", name+"_minus")
        cmd9 = """bdg2bw %s /scratch/sh8tv/Data/Genome/hg38/hg38_clean.len"""%(name+"_plus_FxR.bdg")
        cmd10 = """bdg2bw %s /scratch/sh8tv/Data/Genome/hg38/hg38_clean.len"""%(name+"_minus_FxR.bdg")
        cmd11 = """python /scratch/sh8tv/Script/general/get_signal_gene_Tss_average_NEW.py -w %s,%s,%s,%s -o %s_TSS3kbsig.bed --ext 3000 -i /scratch/sh8tv/Project/regulation_network/Data/Results/hg38_uniq_symbol_withExonArrayExp_TSS.bed"""%(name+"_plus.bw",name+"_minus.bw",name+"_plus_FxR.bw",name+"_minus_FxR.bw",name)
        if readcount >= 10000:
            outf = open('%s.slurm'%(name),'w')
        else:
            outf = open('lowreads/%s.slurm'%(name),'w')
        module = """#!/bin/bash
#SBATCH -n 2
#SBATCH --mem=10GB
#SBATCH -t 18:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o %s.out
#SBATCH -e %s.err
#Run program
"""%(name,name)
        outf.write(module+"\n")
        outf.write(cmd1 + "\n")
        outf.write(cmd2 + "\n")
        outf.write(cmd3 + "\n")
        outf.write(cmd4 + "\n")
        outf.write(cmd5 + "\n")
        outf.write(cmd6 + "\n")
        outf.write(cmd7 + "\n")
        outf.write(cmd8 + "\n")
        outf.write(cmd9 + "\n")
        outf.write(cmd10 + "\n")
        outf.write(cmd11 + "\n")
        outf.close()

