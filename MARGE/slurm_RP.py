import os
allf = os.listdir('../ENCODE_DNase/')
for f in allf:
    if os.path.isfile('../ENCODE_DNase/'+f) and f.endswith('.bw') :#and not f.endswith('_raw.bam'):
        ID = f.split(".")[0]
        model = """#!/bin/bash
#SBATCH -n 8
#SBATCH --mem=100000
#SBATCH -t 10:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o %s.out
#SBATCH -e %s.err
#Run program

python RPCalc.py -b ../ENCODE_DNase/%s.bw -n %s_RP.bed

"""%(ID,ID,ID,ID)
        outf = open('%s.slurm'%(ID),'w')
        outf.write(model)
        outf.close()
        os.system('sbatch %s.slurm'%(ID))
