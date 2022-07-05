import os
allf = os.listdir('.')
for f in allf:
    if os.path.isfile(f) and f.endswith('.bam') and not f.endswith('_raw.bam'):
        ID = f.split(".")[0]
        model = """#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 30:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o %s.log
#SBATCH -e %s.log
#Run program

/nv/vol190/zanglab/sh8tv/Project/scATAC/Data/nakedDNA/Yeast_naked_SRA054922/./downmapPE %s sacCer3
"""%(ID,ID,ID)
        outf = open('%s.slurm'%(ID),'w')
        outf.write(model)
        outf.close()
#        os.system('sbatch %s.slurm'%(ID))
