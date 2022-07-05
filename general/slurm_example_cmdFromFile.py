#import os
inf = open('cmd.sh')
for line in inf:
    SRR = line.split()[1]
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
"""%(SRR,SRR,SRR,'hg38')
    outf=  open('%s.slurm'%(SRR),'w')
    outf.write(model)
    outf.close()
#    os.system("sbatch %s.slurm"%(SRR))

