fqnames=["CAT6_Dfre","CAT6_Dfro","CAT7-d-D","CAT7-f-D","CAT8-f-AT","CAT8-f-D"]

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

/usr/bin/python2.7 /scratch/sh8tv/Script/single_cell/Tanglab/Tang_scATAC_splitBC_pipeline_gzipVer.py -m /nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/meta/meta_coassay.txt -i %s -f /nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/fastq/070809_CAT678/combined_fastq/

"""%(fqname,fqname,fqname)
    outf=  open('%s.slurm'%(fqname),'w')
    outf.write(model)
    outf.close()
    os.system('sbatch %s.slurm'%(fqname))
