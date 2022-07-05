import os,sys
def make_nmer_dict(n):
    nmer_seq = {}
    bp = ['A','C','G','T']
    allseq = [0]*n
    allseq[0] = bp
    i=1
    while i < n:
        allseq[i] = []
        for previous_seq in allseq[i-1]:
            for add_bp in bp:
                new_seq = previous_seq + add_bp
                allseq[i].append(new_seq)
        i += 1
    for seq in allseq[n-1]:
  #      print seq
        nmer_seq[seq] = []
    del allseq
    return nmer_seq

def addbias(encodingdict, infile):
#    os.system('wc -l %s'%infile)
    inf = open(infile)
    tmpdict = {}
    for line in inf:
        ll = line.split()
        seq = ll[0]
        raw = ll[1]
        encoding = ll[2]
        tmpdict[seq] = ll[1:]
    inf.close()
    for seqtype in encodingdict.keys():
        if tmpdict.has_key(seqtype):
            #rawdict[seqtype].append( str(round(float(tmpdict[seqtype][0]),4)) )
            encodingdict[seqtype].append(str(round(float(tmpdict[seqtype][1]),4)) )
        else:
            #rawdict[seqtype].append(str(-8))
            encodingdict[seqtype].append(str(-8))
    return
    
            

### peak bias
Encodingbias = make_nmer_dict(8)
total_sample_list = []

inf = open('biasMat_ownpeak.txt')
outf = open('peak_enc8mer.txt','w')
for line in inf:
    if line.startswith('file'):
        continue
    ll = line.split()
    colname = ll[1]#"_".join(ll[1].split("_")[1:])
    filename = ll[0]
    if os.path.isfile(filename):
        addbias(Encodingbias, filename )
        total_sample_list.append(colname)
    else:
        print 'no file: '+filename
inf.close()

outf.write("\t".join(['seqtype']+total_sample_list)+"\n")
for SeqT in sorted(Encodingbias.keys()):
    newll = [SeqT] + Encodingbias[SeqT]
    outf.write("\t".join(newll)+"\n")
outf.close()

### outpeak bias
Encodingbias = make_nmer_dict(8)
total_sample_list = []

inf = open('biasMat_outpeak.txt')
outf = open('outpeak_enc8mer.txt','w')
for line in inf:
    if line.startswith('file'):
        continue
    ll = line.split()
    colname = ll[1]#"_".join(ll[1].split("_")[1:])
    filename = ll[0]
    if os.path.isfile(filename):
        addbias(Encodingbias, filename )
        total_sample_list.append(colname)
    else:
        print 'no file: '+filename
inf.close()

outf.write("\t".join(['seqtype']+total_sample_list)+"\n")
for SeqT in sorted(Encodingbias.keys()):
    newll = [SeqT] + Encodingbias[SeqT]
    outf.write("\t".join(newll)+"\n")
outf.close()

### chrM bias
Encodingbias = make_nmer_dict(8)
total_sample_list = []

inf = open('biasMat_chrM.txt')
outf = open('chrM_enc8mer.txt','w')
for line in inf:
    if line.startswith('file'):
        continue
    ll = line.split()
    colname = ll[1]#"_".join(ll[1].split("_")[1:])
    filename = ll[0]
    if os.path.isfile(filename):
        addbias(Encodingbias, filename )
        total_sample_list.append(colname)
    else:
        print 'no file: '+filename
inf.close()

outf.write("\t".join(['seqtype']+total_sample_list)+"\n")
for SeqT in sorted(Encodingbias.keys()):
    newll = [SeqT] + Encodingbias[SeqT]
    outf.write("\t".join(newll)+"\n")
outf.close()

### promoter bias
Encodingbias = make_nmer_dict(8)
total_sample_list = []

inf = open('biasMat_promoter3kb.txt')
outf = open('promoter3kb_enc8mer.txt','w')
for line in inf:
    if line.startswith('file'):
        continue
    ll = line.split()
    colname = ll[1]#"_".join(ll[1].split("_")[1:])
    filename = ll[0]
    if os.path.isfile(filename):
        addbias(Encodingbias, filename )
        total_sample_list.append(colname)
    else:
        print 'no file: '+filename
inf.close()

outf.write("\t".join(['seqtype']+total_sample_list)+"\n")
for SeqT in sorted(Encodingbias.keys()):
    newll = [SeqT] + Encodingbias[SeqT]
    outf.write("\t".join(newll)+"\n")
outf.close()

### naked bias
Encodingbias = make_nmer_dict(8)
total_sample_list = []

inf = open('biasMat_naked.txt')
outf = open('naked_enc8mer.txt','w')
for line in inf:
    if line.startswith('file'):
        continue
    ll = line.split()
    colname = ll[1]#"_".join(ll[1].split("_")[1:])
    filename = ll[0]
    if os.path.isfile(filename):
        addbias(Encodingbias, filename )
        total_sample_list.append(colname)
    else:
        print 'no file: '+filename
inf.close()

outf.write("\t".join(['seqtype']+total_sample_list)+"\n")
for SeqT in sorted(Encodingbias.keys()):
    newll = [SeqT] + Encodingbias[SeqT]
    outf.write("\t".join(newll)+"\n")
outf.close()

### genome bias

Encodingbias = make_nmer_dict(8)
total_sample_list = []
folder = '/nv/vol190/zanglab/sh8tv/Project/scATAC/Data/bias_matrix/cutCount/genome/'

inf = open('biasMat_chrM.txt')
outf = open('genome_enc8mer.txt','w')
for line in inf:
    if line.startswith('file'):
        continue
    ll = line.split()
    colname = 'genome_'+"_".join(ll[1].split("_")[1:])
    filename = folder+colname+'_Enc8mer.txt'
    if os.path.isfile(filename):
        addbias(Encodingbias, filename )
        total_sample_list.append(colname)
    else:
        print 'no file: '+filename
inf.close()

outf.write("\t".join(['seqtype']+total_sample_list)+"\n")
for SeqT in sorted(Encodingbias.keys()):
    newll = [SeqT] + Encodingbias[SeqT]
    outf.write("\t".join(newll)+"\n")
outf.close()




