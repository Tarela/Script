fqnames=["CAT6_Rfro","CAT6_Rfre","CAT7-d-R","CAT7-f-R","CAT8-R"]

import os
for fqname in fqnames:
    #fqname = fqnames[i]
    model="""#!/bin/bash
#SBATCH -n 2
#SBATCH --mem=30GB
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o %s.log
#SBATCH -e %s.log
#Run program

module load cutadapt
cutadapt -u 30 -a "A{8}" -n 4 -g AAGCAGTGGTATCAACGCAGAGTACAT -o /nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/fastq/070809_CAT678/combined_fastq/%s_R1trim.fastq.gz /nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/fastq/070809_CAT678/combined_fastq/%s_R1.fastq.gz
/usr/bin/python2.7 /scratch/sh8tv/Script/single_cell/Tanglab/Tang_scRNA_splitBC_pipeline.py -m /nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/meta/meta_coassay.txt -i %s -f /nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/fastq/070809_CAT678/combined_fastq/

"""%(fqname,fqname,fqname,fqname,fqname)
    outf=  open('%s.slurm'%(fqname),'w')
    outf.write(model)
    outf.close()
    os.system('sbatch %s.slurm'%(fqname))

