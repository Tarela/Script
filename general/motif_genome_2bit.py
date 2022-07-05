import sys,os,unittest
import numpy
import seq.count
import seq._seq
import cistrome.regions as cis
import twobitreader

## scan motif whole genome , need binoch package, numpy , input like:
#### python motif_genome.py hg19 hg19.len 20000 peak.bed pwm.txt output_motifsite.bed
def test(seqx,xpssm,out1):#seqs,pssm):
        [xchrom,xstart,xend,xseq] = seqx
        #seqs = ["ACTGATATACTGATACTGCGCGCCAGTGCGCGCGACTGACCAGTACACACACACACC"]
        #["ATATATATATACGCGCGCGCGCGCGACACACACACACACACACACACAC",\
        #"ATATATATATACTCGCGCGCGCGCGACACACACACACACACACACACAC",\
        #"ATATATATATACTGGCGCGCGCGCGACACACACACACACACACACACAC",\
        #"ACTGGGGTACTGAGGCGCGCGCGCGACACACACACACACACACACACAC",\
        #"ATATATATATACGCGCGCGCGCGCGACACACAGTACACACACCCCCAGT",\
        #"ACTGATATACTGATACTGCGCGCGCGCGCGACACACACACACACACACC",\
        #"ACTGATATACTGATACTGCGCGCGCGCGCGACTGACACACACACACACC"]
        #prob = seq.count.count(seqs)
        option = seq._seq.CUTOFF_OPTION
        order = 2
        pssm = numpy.array( xpssm ,numpy.float )
        outf1 = open(out1,'w')
       # outf2 = open(out2,'w')
        #start_t    = [0]#[26,10,10,0,45,0,0]
        #end_t      = [4]#[30,14,14,4,49,4,4]
        #strand_t   = [0]#[1,0,0,0,1,0,0]
        #pssm = numpy.array( [ [1,0,0,0], [0,1,0,0], [0,0,0,1], [0,0,1,0] ] ,numpy.float )i
        for i in range(len(xstart)):
            #print "new 10000 start"
            chrm = xchrom[i]
            start = xstart[i]
            end = xend[i]
            #strand = seqx.strand[i]
            seqs = [xseq[i]]
            #print seqs
            idxmax,startmax,endmax,orientmax,scoremax = seq._seq.seqscan( seqs, pssm, [], order, option )
            #print seq._seq.seqscan( seqs, pssm, [], order, option )
            #exit(1)#print strand
            for j in range(len(startmax)):
                if startmax[j] < 50 :
                    continue 
                realstart = start + startmax[j]
                realend = start + endmax[j]
                if end - realend < 50:
                    continue
                realscore = str(scoremax[j])
                realorient = orientmax[j]
                realseq = seqs[0][int(startmax[j]):int(endmax[j])]
                #ll = [chrm,str(realstart),str(realend),realseq,realscore]
            #    print realseq
            #    print realorient
            #exit(1)
             #   print ll
# for i in range(len(startmax)):
        #    print [str(startmax[i]),str(endmax[i]),str(scoremax[i]),str(orientmax[i]),str(idxmax[i])]
                #newline = "\t".join(ll)+"\n"
                if float(realscore) < 1000:
                    continue
                if realorient == 0:
                    
                    ll = [chrm,str(realstart),str(realend),realseq,realscore,"+"]
                    newline = "\t".join(ll)+"\n"
                    outf1.write(newline)
                elif realorient == 1:
                    
                    ll = [chrm,str(realstart),str(realend),realseq,realscore,"-"]
                    newline = "\t".join(ll)+"\n"
                    outf1.write(newline)
                else:
                    print "error"
           # exit(1)
        outf1.close()
     #   outf2.close()
def makeseq(chromlen,binlength,seq2bit):
    sequence  = seq2bit
    genome = twobitreader.TwoBitFile(sequence) 

    xchrom = []
    xstart = []
    xend = []
    xseq = []
#    x.strand = []
    inf = open(chromlen)#open("/mnt/data/static_libraries/chromLen/hg19.len")
    #outf = open("hg19_%s_bin.bed"%(binlength),'w')
    for line in inf:
        if line.split()[0] == "chrM" or "_" in line.split()[0]:
            continue
        start = 1
        end = int(binlength)
        while end < int(line.split()[1]):
        #    print start,end
            #outf.write("\t".join([line.split()[0],str(start),str(end)])+"\n")
            xchrom.append(line.split()[0])
            xstart.append(max(1,start-100))
            #print line
            xend.append(end+100)
            #print line.split()[0],max(1,start-100),(end+100)
            xseq.append(genome[line.split()[0]][max(1,start-100):(end+100)])
#            x.end.append(min(end+100,int(line.split()[1])))
            start += int(binlength)
            end += int(binlength)
        xchrom.append(line.split()[0])
        xstart.append(max(1,start-100))
        xend.append(int(line.split()[1]))
       # print line.split()[0],max(1,start-100),int(line.split()[1])
        xseq.append(genome[line.split()[0]][max(1,start-100):int(line.split()[1])])

    x = [xchrom,xstart,xend,xseq]
    return x
    #print len( x.seq)
    #outf.write("\t".join([line.split()[0],str(start),line.split()[1]])+"\n")
def makepssm(pssmfile):
    pssm=[]
    inf = open(pssmfile)
    for line in inf:
        if line.startswith("A"):
            continue
        a=line.split()
        if len(a)!=4:
            print a
        else:
            pssm.append(a)
    return pssm
print sys.argv
x=makeseq(sys.argv[2],sys.argv[3],sys.argv[1])
pssm= makepssm(sys.argv[4])
test(x,pssm,sys.argv[5])
