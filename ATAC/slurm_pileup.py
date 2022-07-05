import os
allf = os.listdir('.')
for f in allf:
    if os.path.isfile(f) and f.endswith('.bed') :
        ID = f.split(".")[0]
        model = """#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20GB
#SBATCH -t 20:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o %s.out
#SBATCH -e %s.err
#Run program

#macs2 callpeak -t %s.bed -n %s --nomodel --extsize 50 --keep-dup 1 -q 0.01 -B --SPMR -g hs
#mv %s_treat_pileup.bdg %s.bdg
macs2 pileup -i %s.bed -o %s.bdg -f BED --extsize 1
bdg2bw %s.bdg /home/sh8tv/Main/Data/Genome/hg38/hg38_clean.len

"""%(ID,ID,ID,ID,ID,ID,ID,ID,ID)
        outf = open('%s_pileup.slurm'%(ID),'w')
        outf.write(model)
        outf.close()
        os.system('sbatch %s_pileup.slurm'%(ID))
