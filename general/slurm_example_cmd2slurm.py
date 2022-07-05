import os

def cmd2slurm(CMD):
    ID = CMD.split("-w")[0].split("-o")[-1].strip()[:-4]
    model = """#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=1GB
#SBATCH -t 10:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o %s.log
#SBATCH -e %s.log
#Run program

%s

"""%(ID,ID,CMD)
    OUTF = open('%s.slurm'%(ID),'w')
    OUTF.write(model)
    OUTF.close()
    os.system('sbatch %s.slurm'%(ID))

cmd1 = "python /scratch/sh8tv/Script/general/get_signal_fulllen_average_NEW.py -i /scratch/sh8tv/Project/scATAC/Data/ATAC/human_cellline/GM12878_SEATAC_summit200repov.bed -o GM12878_SEATAC_summit200_GMrep1sig.bed -w SRR891268.bw --bwfolder /nv/vol190/zanglab/sh8tv/Project/scATAC/Data/otherATAC/human_cellline/GM12878SE/"
cmd2 = "python /scratch/sh8tv/Script/general/get_signal_fulllen_average_NEW.py -i /scratch/sh8tv/Project/scATAC/Data/ATAC/human_cellline/GM12878_SEATAC_summit200repov.bed -o GM12878_SEATAC_summit200_GMrep2sig.bed -w SRR891269.bw --bwfolder /nv/vol190/zanglab/sh8tv/Project/scATAC/Data/otherATAC/human_cellline/GM12878SE/"
cmd3 = "python /scratch/sh8tv/Script/general/get_signal_fulllen_average_NEW.py -i /scratch/sh8tv/Project/scATAC/Data/ATAC/human_cellline/GM12878_SEATAC_summit200repov.bed -o GM12878_SEATAC_summit200_K562rep1sig.bed -w SRR2085918.bw --bwfolder /nv/vol190/zanglab/sh8tv/Project/scATAC/Data/otherATAC/human_cellline/K562SE/"
cmd4 = "python /scratch/sh8tv/Script/general/get_signal_fulllen_average_NEW.py -i /scratch/sh8tv/Project/scATAC/Data/ATAC/human_cellline/GM12878_SEATAC_summit200repov.bed -o GM12878_SEATAC_summit200_K562rep2sig.bed -w SRR2085919.bw --bwfolder /nv/vol190/zanglab/sh8tv/Project/scATAC/Data/otherATAC/human_cellline/K562SE/"

cmd2slurm(cmd1)
cmd2slurm(cmd2)
cmd2slurm(cmd3)
cmd2slurm(cmd4)

