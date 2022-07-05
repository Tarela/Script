import math,numpy,sys
from numpy import linalg as la

# A = numpy.array( [ [1,1,-1,-1], [1,-1,1,-1], [1,-1,-1,1], [1,1,1,1] ] ) 
# B = numpy.kron(A,A)

class encoding(object):
    def __init__(self):
        
        self.mononuc = ['A','C','G','T'] 
        self.dinuc   = ['AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT']
 
        self.d1  = { 'A': [1,-1,-1], 'C': [-1,1,-1], 'G': [-1,-1,1], 'T':[1,1,1] }
        self.d2  = { 
            'AA': [+1,-1,-1,-1,+1,+1,-1,+1,+1], 
            'AC': [-1,+1,-1,+1,-1,+1,+1,-1,+1],
            'AG': [-1,-1,+1,+1,+1,-1,+1,+1,-1],
            'AT': [+1,+1,+1,-1,-1,-1,-1,-1,-1],
            'CA': [-1,+1,+1,+1,-1,-1,-1,+1,+1],
            'CC': [+1,-1,+1,-1,+1,-1,+1,-1,+1],
            'CG': [+1,+1,-1,-1,-1,+1,+1,+1,-1], 
            'CT': [-1,-1,-1,+1,+1,+1,-1,-1,-1],
            'GA': [-1,+1,+1,-1,+1,+1,+1,-1,-1],
            'GC': [+1,-1,+1,+1,-1,+1,-1,+1,-1],
            'GG': [+1,+1,-1,+1,+1,-1,-1,-1,+1], 
            'GT': [-1,-1,-1,-1,-1,-1,+1,+1,+1],
            'TA': [+1,-1,-1,+1,-1,-1,+1,-1,-1],
            'TC': [-1,+1,-1,-1,+1,-1,-1,+1,-1],
            'TG': [-1,-1,+1,-1,-1,+1,-1,-1,+1], 
            'TT': [+1,+1,+1,+1,+1,+1,+1,+1,+1] }
            #'GC': [+1,-1,+1,+1,-1,+1,+1,-1,+1],
            #'TG': [-1,-1,+1,-1,-1,+1,-1,-1,-1], 

        # parameters
        self.b0 = 0
        self.b1 = []
        self.b2 = []

        self.B1 = numpy.array( [
            [+1,-1,-1],
            [+1,-1,-1],
            [+1,-1,-1],
            [+1,-1,-1],
            [-1,+1,-1], 
            [-1,+1,-1],
            [-1,+1,-1], 
            [-1,+1,-1], 
            [-1,-1,+1], 
            [-1,-1,+1], 
            [-1,-1,+1], 
            [-1,-1,+1], 
            [+1,+1,+1], 
            [+1,+1,+1], 
            [+1,+1,+1],
            [+1,+1,+1]] 
        )

        self.B2 = numpy.array([
            [+1,-1,-1],
            [-1,+1,-1],
            [-1,-1,+1],
            [+1,+1,+1], 
            [+1,-1,-1],
            [-1,+1,-1],
            [-1,-1,+1],
            [+1,+1,+1], 
            [+1,-1,-1],
            [-1,+1,-1],
            [-1,-1,+1],
            [+1,+1,+1], 
            [+1,-1,-1],
            [-1,+1,-1],
            [-1,-1,+1],
            [+1,+1,+1] ]
        )
 
        self.H = numpy.array( 
            [[+1,-1,-1,-1,+1,+1,-1,+1,+1],
             [-1,+1,-1,+1,-1,+1,+1,-1,+1],
             [-1,-1,+1,+1,+1,-1,+1,+1,-1],
             [+1,+1,+1,-1,-1,-1,-1,-1,-1],
             [-1,+1,+1,+1,-1,-1,-1,+1,+1],
             [+1,-1,+1,-1,+1,-1,+1,-1,+1],
             [+1,+1,-1,-1,-1,+1,+1,+1,-1],
             [-1,-1,-1,+1,+1,+1,-1,-1,-1],
             [-1,+1,+1,-1,+1,+1,+1,-1,-1],
             [+1,-1,+1,+1,-1,+1,-1,+1,-1],
             [+1,+1,-1,+1,+1,-1,-1,-1,+1],
             [-1,-1,-1,-1,-1,-1,+1,+1,+1],
             [+1,-1,-1,+1,-1,-1,+1,-1,-1],
             [-1,+1,-1,-1,+1,-1,-1,+1,-1],
             [-1,-1,+1,-1,-1,+1,-1,-1,+1],
             [+1,+1,+1,+1,+1,+1,+1,+1,+1]] )

    def check(self):
        A = numpy.hstack( (self.B1, numpy.hstack( (self.B2, self.H) ) ))
        #print numpy.dot( numpy.transpose( A ), A ) 

    def rparam(self):
        """
        compute parameters for reverse orientation
        """
        rcmono = {'A':'T','C':'G','G':'C','T':'A'} 
        rcdi   = {'AA':'TT', 'AC':'GT', 'AG':'CT', 'AT':'AT', 'CA':'TG', 'CC':'GG', 'CG':'CG', 'CT':'AG', 'GA':'TC', 'GC':'GC', 'GG':'CC', 'GT':'AC', 'TA':'TA', 'TC':'GA', 'TG':'CA', 'TT':'AA'}
         
        A = numpy.hstack( ( numpy.transpose(numpy.array([16*[1]])), numpy.hstack( (self.B1, numpy.hstack( (self.B2, self.H) ) ))) )
        #print numpy.dot( numpy.transpose( A ), A ) 

        R = []
        for i,elem in enumerate(self.dinuc):
            rc = rcdi[elem]
            s = [1] + self.d1[rc[0]] + self.d1[rc[1]] + self.d2[rc]   
            R.append(s)  
        #print  numpy.array(R) 
        #print numpy.dot( numpy.transpose( R ), R ) 
        self.RCMAP =  numpy.dot( 1.0/16*numpy.transpose( A ), R ) 


def parseseqratio(file):
    fp = open(file)
    seqd = {}
    seqr = {}
    for line in fp.readlines():
        f = line.split()
        try:
            seqd[f[0]] = math.log(float(f[1]))
            seqr[f[0]] = math.log(float(f[2]))
        except:
            pass
    return seqd,seqr


def paramest(seqd):
    """
    Implementation of
    Maximally Efficient Modeling of DNA Sequence Motifs at All Levels of Complexity
    Gary D. Stormo
    Genetics, 2011
    """

    code=encoding()
    #code.check()
    #code.rparam()

    mat,y  = [],[]

    for j,seq in enumerate(seqd):
        l = list(seq)
    
        # constant
        s = [1]
        
        # single nucleotides 
        for i,elem in enumerate( l ):
            s += code.d1[l[i]] 
      
        # interaction of adjacent nucleotides 
        for i,elem in enumerate( l[0:-1] ):
            s += code.d2[ l[i]+l[i+1] ]
 
        mat.append(s)
        y.append(seqd[seq])
    y = numpy.array(y)
    b = ols( numpy.array(mat), y )
    #print b
    b0 = b[0]
    b1 = b[1:(1+3*len(seq))]
    b1 = numpy.reshape(b1,(-1,3))
    b2 = b[(1+3*len(seq)):len(b)]
    b2 = numpy.reshape(b2,(-1,9))
    return b0,b1,b2


def ols(X,y):
    # b=(X'X)^-1X'y
    XT = numpy.transpose(X)
    Z = numpy.dot(XT,X) 
    A = numpy.dot( la.inv( numpy.dot(XT,X) ), XT )
    b = numpy.dot(A,y)

    yp =  numpy.dot(X,b)
    #for i,elem in enumerate(y):
    #   print y[i],yp[i]
    return b


def predict2(seq,b0,b1,b2,w=6):
    """
    prediction based on w-mer using 
    mono- and di-nucleotides
    """
    n = len(seq)
    code = encoding()
    mono = numpy.array( list(seq) )
    di   = numpy.array( [ seq[i]+seq[i+1] for i,elem in enumerate(seq[0:-1]) ] ) 

    # make design matrix from sequence
    # mono-
    X1 = numpy.zeros( (len(mono),3) )
    # di-
    X2 = numpy.zeros( (len(di),9)   )

    for nuc,v in code.d1.iteritems():
        X1[mono==nuc] = numpy.array(v)
    for dinuc,v in code.d2.iteritems(): 
        X2[di==dinuc] = numpy.array(v)

    Y0 = b0*numpy.ones(len(mono))
    Y1 = numpy.dot(X1,numpy.transpose(b1))
    Y2 = numpy.dot(X2,numpy.transpose(b2))

    Y = Y0[0:(n-(w-1))]
    for i in range(w):
        Y += Y1[i:(n-(w-1)+i),i] 

    for i in range(w-1):
        Y += Y2[i:(n-(w-1)+i),i]

    return Y


def revcomp(seq):
    rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A' }
    seqt = list(seq)
    seqt.reverse()
    r = ''.join( [ rc[x] for x in seqt] )
    return r


def ATAC(shift0file,shift2file,seqlist,offset=8):
    code = encoding() 
    seqd,seqr = parseseqratio(shift0file)
    b0s0,b1s0,b2s0 = paramest(seqd)

    #for elem in seqd.keys():
    #    yf = predict2(elem,b0s0,b1s0,b2s0,w=6)
    #    print elem,seqd[elem],yf[0] 

    seqd,seqr = parseseqratio(shift2file)
    #print seqd
    b0s2,b1s2,b2s2 = paramest(seqd)

    #for elem in seqd.keys():
    #    yf = predict2(elem,b0s2,b1s2,b2s2,w=6)
    #    print elem,seqd[elem],yf[0] 

    # forward
    b0f = 6.0/8*b0s0 + 2.0/8*b0s2 

    # mononucleotides - concatenate parameters from regression on 6-mers 
    b1f = numpy.vstack( (b1s0, b1s2[3:-1,:] ) )
    # last nucleotide is counted twice so multiply by 0.5
    b1f = numpy.vstack( (b1f, 0.5*b1s2[-1,:]) )
    # concatenate dinucleotides
    b2f = numpy.vstack( (b2s0,b2s2[3:,:])     )
    print b0f,b1f,b2f
    
    for seq in seqlist: 
        # predict forward strand
        yf  = predict2(seq,b0f,b1f,b2f,w=8)
 
        # predict reverse strand
        yr  = predict2(revcomp(seq),b0f,b1f,b2f,w=8)
    
        # reverse reverse prediction
        yrr = yr[::-1] 
        #print yr
        #print yrr
        # add forward and reverse with offset
        z   = yf[:-offset] + yrr[offset:]

        #for elem in z:
        #    print elem
    return

 
if __name__ == "__main__":

    # 6 mer ratio file, 0 shift
    shift0file = 'GM12878_ATACseq_50k_bias_onPeak_u3d3s0_newcode.txt'

    # 6 mer ratio file, 2 shift
    shift2file = 'GM12878_ATACseq_50k_bias_onPeak_u3d3s2_newcode.txt'
    seqlist    = ['ATATATATCACCTCCAGGGGACACACCAATTTTTGATGATTTGCGCTGT','GCGCGCTATATTATAGGCGCGTGATCGTCACGCGCGCGAGAGAGTTTTATATGGGCG']
    ATAC(shift0file,shift2file,seqlist)

