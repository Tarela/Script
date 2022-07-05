import sys,os,unittest
import numpy
import seq.count
import seq._seq
import cistrome.regions as cis
### test5 is on given region
def test(seqx,xpssm,out1):#seqs,pssm):

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
        for i in range(len(seqx.start)):
            #print "new 10000 start"
            chrm = seqx.chrom[i]
            start = seqx.start[i]
            #strand = seqx.strand[i]
            seqs = [seqx.seq[i]]
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
                if float(realscore) < 1000 : 
                    continue
                if realorient == 0:
                    ll = [chrm,str(realstart-1),str(realend-1),realseq,realscore,"+"]
                    newline = "\t".join(ll)+"\n"
                    outf1.write(newline)
                elif realorient == 1:
                    ll = [chrm,str(realstart-2),str(realend-2),realseq,realscore,"-"]
                    newline = "\t".join(ll)+"\n"
                    outf1.write(newline)
                else:
                    print "error"
           # exit(1)
        outf1.close()
     #   outf2.close()
def makeseq(bedfile,species):
    x = cis.interval(genome=species)
    x.chrom = []
    x.start = []
    x.end = []
#    x.strand = []
    inf = open(bedfile)#open("/mnt/data/static_libraries/chromLen/hg19.len")
    #outf = open("hg19_%s_bin.bed"%(binlength),'w')
    for line in inf:
        #start = 1
        #end = int(binlength)
        #while end < int(line.split()[1]):
        #    print start,end
            #outf.write("\t".join([line.split()[0],str(start),str(end)])+"\n")
        x.chrom.append(line.split()[0])
        x.start.append(int(line.split()[1]))
        x.end.append(int(line.split()[2]))
        #    start += int(binlength)
        #    end += int(binlength)
        #x.chrom.append(line.split()[0])
        #x.start.append(start)
        #x.end.append(int(line.split()[1]))
    x.getSequence()
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
x=makeseq(sys.argv[2],sys.argv[1])
pssm= makepssm(sys.argv[3])
test(x,pssm,sys.argv[4])
