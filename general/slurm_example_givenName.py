SRRs=['SRR1533850','SRR1533851','SRR1533852','SRR2085919','SRR2096437']
SPEs=['mm10','mm10','mm10','hg38','hg38']
import os
for i in range(5):
    SRR = SRRs[i]
    species = SPEs[i]
    model="""#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 10:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o %s.log
#SBATCH -e %s.log
#Run program

./downmap %s %s
"""%(SRR,SRR,SRR,species)
    outf=  open('%s.slurm'%(SRR),'w')
    outf.write(model)
    outf.close()
    os.system('sbatch %s.slurm'%(SRR))
