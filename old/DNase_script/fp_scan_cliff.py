import argparse,sys,numpy
from cistrome import regions as c

PLUS       = 0
MINUS      = 1
MOTIFCUT   = 1e4 
#MOTIFCUT  = 0
DEBUG      = True
MISSINGVAL = -99


def compute_DHS_motif( motiflen, cuts_positive, cuts_negative ):
    """
    Compute motif from cuts.
    """
    cp = {}
    cn = {}
    #print cuts_positive
    #print cuts_positive['+'].sum(0)
    #print cuts_positive['-'].sum(0)
    #print cuts_negative['+'].sum(0)
    #print cuts_negative['-'].sum(0)
    #print motiflen
    cp['+'] = cuts_positive['+'].sum(0)
    #cp['-'] = cuts_positive['-'].sum(0)
    cn['+'] = cuts_negative['+'].sum(0)
    #cn['-'] = cuts_negative['-'].sum(0)

    tot = cp['+'].sum() + cn['+'].sum() 
    DHS_motif = ( cp['+']/tot, cn['+']/tot )
    print 'motif+:', ','.join( [ '%5.3f' % x for x in DHS_motif[0] ] )
    print 'motif-:', ','.join( [ '%5.3f' % x for x in DHS_motif[1] ] )
    return DHS_motif


def compute_DHS_modelscore( cuts_positive, cuts_negative, DHSval_p, DHSval_n, span=20, pseudocount=0.5 ):
    """
    Compute motif from cuts.
    """

    EPS = 1e-3

    # check these are all the same length
    assert( len( cuts_positive ) == len( cuts_negative ) )
    assert( len( cuts_positive ) == len( DHSval_p ) )
    assert( len( cuts_positive ) == len( DHSval_n ) )

    n = len(cuts_positive)

    DHSval_n = numpy.array(DHSval_n)
    DHSval_p = numpy.array(DHSval_p)

    score   = []
    score_n = []
    score_p = []
    tot_n   = []
    tot_p   = []

    # loop over length of sequence
    for i in range( n - span ):

        observe_plus  = numpy.array( cuts_positive[i:(i+span)] )
        observe_minus = numpy.array( cuts_negative[i:(i+span)] )
        observe_sum   = sum(observe_plus) + sum(observe_minus)

        if observe_sum > 0:

            # predicted cuts
            model_plus  = DHSval_p[i:(i+span)]
            model_minus = DHSval_n[i:(i+span)]

            if ( min(model_plus) == MISSINGVAL ) or ( min(model_minus) == MISSINGVAL ):
                f_p    = MISSINGVAL   
                f_n    = MISSINGVAL 
                f_tot  = MISSINGVAL
            else:
                model_plus  += EPS
                model_minus += EPS

                model_sum = sum(model_plus) + sum(model_minus)
                # model_plus    /= model_sum # normalize 
                # model_minus   /= model_sum 
                # treat plus and minus separately
                model_plus  /= sum( model_plus  ) # normalize 
                model_minus /= sum( model_minus )

                # observed cuts with predicted prior 
                lambda_plus  = observe_plus  + pseudocount * model_plus
                lambda_minus = observe_minus + pseudocount * model_minus
                #observe_sum   = sum(observe_plus) + sum(observe_minus) + 2*pseudocount
                observe_sum  = sum(observe_plus) + sum(observe_minus)

                # normalize 
                #lambda_plus   /= observe_sum 
                #lambda_minus  /= observe_sum
                # normalize 
                lambda_plus   /= sum( lambda_plus  ) 
                lambda_minus  /= sum( lambda_minus )

                lambda_plus0  = model_plus
                lambda_minus0 = model_minus

                f_plus  = numpy.sum( observe_plus  * ( numpy.log(lambda_plus)  - numpy.log(lambda_plus0)  ) ) 
                f_minus = numpy.sum( observe_minus * ( numpy.log(lambda_minus) - numpy.log(lambda_minus0) ) )
                #f_tot   = ( f_plus + f_minus )/( observe_sum )
                f_p     = ( f_plus  )/( max(1,sum(observe_plus))  )
                f_n     = ( f_minus )/( max(1,sum(observe_minus)) )
                f_tot   = ( f_p + f_n )
                # print location:
                #print  '\t'.join( ['%s' % loc['+'][i][0], '%d' % loc['+'][i][1], '%d' % loc['+'][i][2], '%s' % loc['+'][i][3], '%5.3f' % loc['+'][i][4], '%5.3f' % f1 ] )+ '\t' + '\t'.join( [ '%d' % x for x in observe_plus] + [ '%d' % x for x in observe_minus] )
            ####print '%4.2f' % f1 #'\t'.join( [ '%4.2f' % x for x in model_plus ] + [ '%4.2f' % x for x in observe_plus ] )
            
            score_n.append(f_n)
            score_p.append(f_p)
            score.append(f_tot)
            tot_n.append( sum(observe_minus) )
            tot_p.append( sum(observe_plus) )

        else:
            f_tot = 0
            score_n.append(0)
            score_p.append(0)
            score.append(0)
            tot_n.append( sum(observe_minus) )
            tot_p.append( sum(observe_plus) )

    score   = numpy.array(score)
    score_n = numpy.array(score_n)
    score_p = numpy.array(score_p)
    imax = numpy.where( score == max(score) )[0]

    return imax[0], score[imax[0]], score, score_n, score_p, tot_n, tot_p 



def read_file(fp,mintags=0,maxtags=10,maxlines=100,select=100):
    """
    Parse data
    """

    # TODO read in more general format
    if 1:
        CHR, START, END, NAME, SCORE, STRAND, COUNT50P, COUNT50M = 0,1,2,3,4,5,6,7
        CUT_START = 14
        CUT_END   = 414
    else:
        CHR, START, END, NAME, SCORE, COUNT50P, COUNT50M = 0,1,2,3,4,5,6
        STRAND = None
        CUT_START = 13
        CUT_END   = 413

    k = 0
    X = c.interval(genome='hg19')
    X.chrom,X.start,X.end,X.strand,X.name,X.val = [],[],[],[],[],[]
    cuts_pos = []
    cuts_neg = []

    mid         = int( 0.5 * ( CUT_END - CUT_START ) )
    startoffset = ( mid - int( 0.5 * select ) )
    endoffset   = ( mid - int( 0.5 * select ) + select )

    for elem in fp.readlines():
        if elem[0:3] == 'chr':
            f = elem.split()
            chr,start,end,seq,motifscore,n50p,n50m = f[CHR], int(f[START]), int(f[END]), f[NAME], float(f[SCORE]), float(f[COUNT50P]), float(f[COUNT50M] )

            if STRAND:
                strand = f[STRAND]
            else: 
                strand = '+'

            if ( strand == '+' ) and ( n50p + n50m >= mintags ) and ( n50p + n50m < maxtags ):
                ## Shawn : only + motif included
                X.chrom.append( f[0] )
                #X.start.append( start - (CUT_END-CUT_START)/2 )
                #X.end.append( start + (CUT_END-CUT_START)/2 )
                X.start.append( start - (CUT_END-CUT_START)/2 + startoffset )
                X.end.append( start - (CUT_END-CUT_START)/2 + endoffset )

                X.strand.append( strand )
                X.name.append( seq )
                X.val.append( (motifscore,n50p+n50m) )
                k += 1

                poscut = [ float(z) for z in f[ CUT_START: CUT_END] ]
                negcut = [ float(z) for z in f[ CUT_END: 2*CUT_END+2-CUT_START ] ]   
                ##Shawn : negcut = [ float(z) for z in f[ CUT_END: 2*CUT_END-CUT_START ] ]

                cuts_pos.append( poscut[ startoffset: endoffset ] ) 
                cuts_neg.append( negcut[ startoffset: endoffset ] )
 

        if k == maxlines:
            break

    X.getSequence()

    #for i,elem in enumerate( X.seq ):
    #    print X.name[i], elem[ (CUT_END-CUT_START)/2: (CUT_END-CUT_START)/2 + 10  ], X.strand[i],  '\t'.join(  [ '%3.1f' % x for x in cuts_pos[i][ (CUT_END-CUT_START)/2: (CUT_END-CUT_START)/2 + 10 ] ] ) 

        #X.seq[i]    = elem[ (CUT_END-CUT_START)/2: (CUT_END-CUT_START)/2 + 10  ]
        #cuts_pos[i] = cuts_pos[i][ (CUT_END-CUT_START)/2: (CUT_END-CUT_START)/2 + 10 ]

    #return
    return X, numpy.array(cuts_pos), numpy.array(cuts_neg)


def count_cut_nmers( seqlist, cuts_pos, cuts_neg, nmer=4, offset=0 ):
    """
    count the number of cuts associated with each nmer in sequence covered by X.
    offset is the position of the cut to be associated with each nmer.
    if offset = 0 the first base of the tag is lined up with the nmer start
    """

    # keep count of the number of occurrences of each n-mer
    allcut_nmer_dict   = {}

    cut_nmer_dict_p    = {}
    nocut_nmer_dict_p  = {}

    cut_nmer_dict_n    = {}
    nocut_nmer_dict_n  = {}

    for i,elem in enumerate(seqlist):

        s  = seqlist[i]
        ns = len(s)
        s  = s.upper()
        ## TODO s = s[0:len(s)/2]
        #s  = s[ (ns/2-flank): (ns/2+flank) ]
        #cp = cuts_pos[ i, (ns/2-flank): (ns/2+flank) ]
        #cn = cuts_neg[ i, (ns/2-flank): (ns/2+flank) ]
        cp = cuts_pos[i,:]
        cn = cuts_neg[i,:]

        #print 'xl', len(X.seq[i]), '\t'.join( [ str(x) for x in cuts_pos[i,0:400] ] ), '\t'.join( [ str(x) for x in cuts_neg[i,0:400] ]  )
        #print 'xs', '\t'.join( [ str(x) for x in cp ] ), '\t'.join( [ str(x) for x in cn ] )

        for k in range( max(0,-offset), len(s) - max(nmer,offset) ):
            #print s[k:k+nmer], cuts_pos[i,k] 

            if 'N' not in list(s[k:k+nmer]):

                if cp[k+offset] > 0:
                    try: 
                        cut_nmer_dict_p[ s[k:k+nmer] ] += cp[k+offset]
                    except:
                        cut_nmer_dict_p[ s[k:k+nmer] ]  = cp[k+offset]
                else:
                    try: 
                        nocut_nmer_dict_p[ s[k:k+nmer] ] += 1
                    except:
                        nocut_nmer_dict_p[ s[k:k+nmer] ]  = 1

                try: 
                    allcut_nmer_dict[ s[k:k+nmer] ] += 1
                except:
                    allcut_nmer_dict[ s[k:k+nmer] ]  = 1

                # negative cuts
                if cn[k+offset-1] > 0:
                    try: 
                        cut_nmer_dict_n[ s[k:k+nmer] ] += cn[k+offset-1]  # note the -1 for the negative strand reads, the cut is to the right
                    except:
                        cut_nmer_dict_n[ s[k:k+nmer] ]  = cn[k+offset-1]  # note: -1
                else:
                    try: 
                        nocut_nmer_dict_n[ s[k:k+nmer] ] += 1
                    except:
                        nocut_nmer_dict_n[ s[k:k+nmer] ]  = 1

    return allcut_nmer_dict, cut_nmer_dict_p, nocut_nmer_dict_p, cut_nmer_dict_n, nocut_nmer_dict_n


def build_nmer_model(cuts_pos,cuts_neg,seqlist,lflank=3,rflank=3):
    """
    count cuts within each nmer sequence 
    """
    # |0.1.2.3.4 lflank=0 rflank=5
    #  0|1.2.3.4 lflank=1 rflank=4
    nmer = lflank + rflank

    # count  
    allcut_nmer_dict, cut_nmer_dict_p, nocut_nmer_dict_p, cut_nmer_dict_n, nocut_nmer_dict_n = count_cut_nmers( seqlist, cuts_pos, cuts_neg, nmer=nmer, offset=lflank )

    dall_p = {}
    dall_n = {}

    # normalize
    for elem in allcut_nmer_dict:

        try:
            val_p = cut_nmer_dict_p[elem]
        except:
            val_p = 0

        try:
            val_n = cut_nmer_dict_n[elem]
        except:
            val_n = 0
 
        try:
            dall_p[elem] += [ val_p / (1.0*allcut_nmer_dict[elem]) ]
            dall_n[elem] += [ val_n / (1.0*allcut_nmer_dict[elem]) ]
        except:
            dall_p[elem]  = [ val_p / (1.0*allcut_nmer_dict[elem]) ]
            dall_n[elem]  = [ val_n / (1.0*allcut_nmer_dict[elem]) ]

    # print
    if 0:
        for elem in dall_p:
            try:
                print elem, allcut_nmer_dict[elem], '\t'.join( [ '%4.3f' % x for x in dall_p[elem] ] ), '\t'.join( [ '%4.3f' % x for x in dall_n[elem] ] )
            except:
                pass

    return dall_p, dall_n


def modelDHS( DHSlookup_p, DHSlookup_n, cuts_positive, cuts_negative, loc, lflank=2, rflank=2, span=20, pseudocount=5 , bdg=None):
    """
    predict DHS cuts and score difference from expectation
    """
    DHSvalP,DHSvalN = [],[]
    seqlist = loc.seq

    for i,seq in enumerate(seqlist):

        DHSval_n,DHSval_p = predictDHSfromSeq( DHSlookup_p, DHSlookup_n, seq, lflank=lflank, rflank=rflank )
        imax,scoremax,score,score_n,score_p,tot_n,tot_p = compute_DHS_modelscore( cuts_positive[i], cuts_negative[i], DHSval_p, DHSval_n, span=span, pseudocount=pseudocount )
        #print '\t'.join( [ '%4.2f' % x for x in score ] + '\t'.join( [ '%4.2f' % x for x in cuts_positive[i] ] + [ '%4.2f' % x for x in cuts_negative[i] ] + [ '%4.2f' % x for x in DHSval_p + DHSval_n ] ) 
        print loc.chrom[i], loc.start[i] + imax, loc.start[i] + imax + span, seq[imax:(imax+span)], loc.val[i][0], loc.val[i][1], scoremax, imax
        #','.join( [ '%d' % x for x in cuts_positive[i][imax:(imax+span)] ] ),\
        #','.join( [ '%d' % x for x in cuts_negative[i][imax:(imax+span)] ] ),\
        #','.join( [ '%4.2f' % x for x in DHSval_p[imax:(imax+span)] ] ),\
        #','.join( [ '%4.2f' % x for x in DHSval_n[imax:(imax+span)] ] )

        if bdg:
            for k,scoret in enumerate(score): 
                print >> bdg, loc.chrom[i], loc.start[i] + k, loc.start[i] + k + span, scoret, score_p[k], score_n[k], tot_n[k], tot_p[k]
    return


def predictDHSfromSeq( DHSlookup_p, DHSlookup_n, seq, lflank=2, rflank=2 ): 
    """
    input sequence string and DHS prediction dictionary
    .0.1.2.3.4.5.
    |0.1.2.3.4.5. lflank=0, rflank=6 
    .0|1.2.3.4.5. lflank=1, rflank=5
    .0.1|2.3.4.5. lflank=2, rflank=4
    .0.1.2|3.4.5. lflank=3, rflank=3  + cut at 3, - cut at 2
    """

    DEFAULTVAL = 0.5

    DHSval_p = []
    DHSval_n = []

    for i in range(len(seq)):
        if ( i >= lflank ) and ( i < len(seq) - rflank ):
            smer_p = seq[(i-lflank)   : (i+rflank)   ]    # + strand 
            smer_n = seq[(i-lflank+1) : (i+rflank+1) ]    # - strand
            try:
                DHSval_p.append( DHSlookup_p[smer_p.upper()][0] )
                DHSval_n.append( DHSlookup_n[smer_n.upper()][0] )
            except: 
                DHSval_p.append( DEFAULTVAL )
                DHSval_n.append( DEFAULTVAL )
        else:
            DHSval_p.append( MISSINGVAL )
            DHSval_n.append( MISSINGVAL )

    return DHSval_n,DHSval_p


def print_nmer_model( DHSlookup_p, DHSlookup_n, fp ):
    #fp = open(fname,'w')
    for elem in DHSlookup_p:
        print >> fp, elem, '%5.3f' % DHSlookup_p[elem][0], '%5.3f' % DHSlookup_n[elem][0]
    #fp.close()
    return

def read_nmer_model( fp ):
    DHSlookup_p,DHSlookup_n = {},{}
    #fp = open(fname)
    for line in fp.readlines():
        f = line.split()
        nmer,pval,nval = f[0],float(f[1]),float(f[2]) 
        DHSlookup_p[nmer] = [pval] 
        DHSlookup_n[nmer] = [nval] 
    fp.close()
    return DHSlookup_p, DHSlookup_n


def test():

    for i in range(10):
        imax, scoremax, score, score_n, score_p, tot_n, tot_p = compute_DHS_modelscore( numpy.array([0,i,0,0]), numpy.array([0,0,0,0]), numpy.array([0.25,0.25,0.25,0.25]), numpy.array([0.25,0.25,0.25,0.25]), span=3, pseudocount=1 )
        print 'imax',imax
        print 'scoremax',scoremax
        print 'score',score
        print 'score_n',score_n
        print 'score_p',score_p
        print 'tot_n',tot_n
        print 'tot_p',tot_p 
        print '========================='
    return 


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
        """Searches for footprints based on deviation from DNA sequence predicted cut biases""")

    parser.add_argument("-i","--input",  dest="input",  default="../SEPE_AR/AR_FPstore.hi.bed", type=argparse.FileType('r'), help="input file" )
    parser.add_argument("-o","--output", dest="output", default=None, type=argparse.FileType('w'), help="output file" )
    parser.add_argument("-c","--count",dest="pseudocount",type=int, default=5,
                         help="Number of pseudocount reads given as prior information (default 5).")
    parser.add_argument("-s","--span",dest="span",type=int,default=20,
                         help="Number of base pairs used in footprint scan (default 20).")
    parser.add_argument("-l","--lflank",dest="lflank",type=int,default=3,
                         help="Number of base pairs on left flank of cut site (default 3).")
    parser.add_argument("-r","--rflank",dest="rflank",type=int,default=3,
                         help="Number of base pairs on right flank of cut site.(default 3).")
    parser.add_argument("--mintags",dest="mintags",type=int,default=50,
                         help="Minimum number of tags in model building region.")
    parser.add_argument("--maxtags",dest="maxtags",type=int,default=500,
                         help="Maximum number of tags in model building region.")
    parser.add_argument("-p","--parameters",dest="paramfile",type=argparse.FileType('r'),
                        help="File of cut bias parameters to be used to predict cut pattern.")
    parser.add_argument("-w","--writeparameters",dest="wparamfile",type=argparse.FileType('w'),
                        help="File of cut bias parameters to record cut bias.")
    parser.add_argument("-x","--xscan",dest="runfpscan",action='store_true',default=False,
                        help="Run foot print scan after generating model")
    parser.add_argument("-m","--maxlines",dest="maxlines",type=int,
                        help="Number of lines to be read from input file.")
    parser.add_argument("--select",dest="select",type=int,
                        help="Number of bases from the center to be analyzed.")

    args = parser.parse_args()

    maxlines    = args.maxlines
    #read_model_from_file = True
    lflank      = args.lflank
    rflank      = args.rflank
    span        = args.span
    pseudocount = args.pseudocount
    mintags     = args.mintags
    maxtags     = args.maxtags

    #option_nmer   = 6
    #option_nmer  = 4
    #fp = open( "../Lncap_SEPE_DHS_50frag_summits_cleavageStore.bed" )   
    #fp = open( "../K562_DHS_50frag_summits_cleavageStore.bed" )
    #fp = open( "../GM06990_DHS_50frag_summits_cleavageStore.bed" )
    #fp = open( "../SEPE_AR/AR_FPstore.hi.bed" )
    fpin  = args.input
    fpout = args.output
    

    #test()
    #sys.exit()

    #fpout = open( "AR_FP.bdg", 'w' )
    X, cuts_pos, cuts_neg = read_file(fpin, mintags=mintags, maxtags=maxtags, maxlines=maxlines, select=args.select )

    # file containing nmer statistics 
    if args.paramfile: 
        DHSlookup_p, DHSlookup_n = read_nmer_model( args.paramfile )
        #DHSlookup_p, DHSlookup_n = read_nmer_model( 'LNCaP_3L_3R.txt' )
    else:
        DHSlookup_p, DHSlookup_n = build_nmer_model(cuts_pos,cuts_neg,X.seq,lflank=lflank,rflank=rflank)
        #print_nmer_model( DHSlookup_p, DHSlookup_n, 'LNCaP_3L_3R.txt' )
        print_nmer_model( DHSlookup_p, DHSlookup_n, args.wparamfile )

    if args.runfpscan == True:
        modelDHS( DHSlookup_p, DHSlookup_n, cuts_pos, cuts_neg, X, lflank=lflank, rflank=rflank, span=span, pseudocount=pseudocount, bdg=fpout )
 
