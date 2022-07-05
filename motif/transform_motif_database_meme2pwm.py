import sys
inf = open(sys.argv[1])
#inf = open('/mnt/Storage/home/huse/Data/Motif/meme/JASPAR_CORE_2016_vertebrates.meme')
#inf = open("/nv/vol190/zanglab/shared/Motif/database/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme")
#inf = open("/nv/vol190/zanglab/shared/Motif/database/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme")
nowWriting = 0
while 1:
    line = inf.readline()
    if line == '':
        break
    if line.startswith('MOTIF'):
 #       print line.strip(), line.split()[-1], line.split()[-1].split("_")[0]
 #       continue
        motifname = line.split()[-1].split("_")[0]
        #motifname = line.split()[2].replace('(','').replace(')','') + "_" + line.split()[1]
        outf = open(motifname+".txt",'w')
        nowWriting = 1
    if line.startswith('letter-probability') and nowWriting ==1:
        #continue
        while 1:
            Lpwm = inf.readline()
         
            if Lpwm.strip() == "" or Lpwm.startswith("URL"):
                nowWriting = 0
                outf.close()
                break
            else:
                newll = Lpwm.strip().split()
                if not len(newll) == 4:
                    print newll,'error motif len'
                outf.write("\t".join(newll) + "\n")
inf.close()


