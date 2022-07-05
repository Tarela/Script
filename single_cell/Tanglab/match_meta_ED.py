import os,sys
from Levenshtein import distance


#edit_dis = distance("ACACAC", "AAAAAT")
P5_refer = ["TATAGCCT","ATAGAGGC"]#,"CCTATCCT","GGCTCTGA","AGGCGAAG","TAATCTTA","CAGGACGT","GTACTGAC"]
P7_refer = ["CGAGTAAT","TCTCCGGA","AATGAGCG","GGAATCTC"]#,"TTCTGAAT","ACGAATTC","AGCTTCAG","GCGCATTA","CATAGCCG","TTCGCGGA","GCGCGAGA","CTATCGCT"]

#P5_refer_use = ["TATAGCCT","ATAGAGGC"]
#P7_refer_use = ["CGAGTAAT","TCTCCGGA","AATGAGCG","GGAATCTC"]

def ED_twopair(B1,B2):
    b1p1 = B1.split("_")[0]
    b1p2 = B1.split("_")[1]
    b2p1 = B2.split("_")[0]
    b2p2 = B2.split("_")[1]
    return "%s:%s"%(distance(b1p1,b2p1),distance(b1p2,b2p2))

def single_barcode_adj(b,BCtype):
    dis_to_refer = []
    if BCtype == 5:
        if b in P5_refer[:2]:
            return b
        else:
            for referBC in P5_refer:
                dis_to_refer.append( distance(b,referBC) )
            sorted_dis_to_refer = sorted(dis_to_refer)
            if sorted_dis_to_refer[0] <= 1 and \
               sorted_dis_to_refer[1] - sorted_dis_to_refer[0] >= 1 and \
               dis_to_refer.index(sorted_dis_to_refer[0]) <= 1:
                return P5_refer[dis_to_refer.index(sorted_dis_to_refer[0])]
            else:
                return "NA"
    if BCtype == 7:
        if b in P7_refer[:4]:
            return b
        else:
            for referBC in P7_refer:
                dis_to_refer.append( distance(b,referBC) )
            sorted_dis_to_refer = sorted(dis_to_refer)
            if sorted_dis_to_refer[0] <= 1 and \
               sorted_dis_to_refer[1] - sorted_dis_to_refer[0] >= 1 and \
               dis_to_refer.index(sorted_dis_to_refer[0]) <= 3:
                return P7_refer[dis_to_refer.index(sorted_dis_to_refer[0])]
            else:
                return "NA"    


cell_barcode_dict = {}
inf = open("meta.txt")
for line in inf:
    ll = line.split()
    cell_barcode_dict[ll[0]+":"+ll[3]] = ll[5]
inf.close()

cell_vs_adjbarcode_dict = {}
for f in os.listdir("./"):
    if f.endswith("_BCcount.txt") and os.path.isfile(f):
        inf = open(f)
        batch_name = f[:-12]
        outf = open("%s_BCex.txt"%batch_name,'w')
        #print batch_name
        if "Dpool" in f:
            topN = 8
        else:
            topN = 1 

        BCnum = 0
        real_barcodes = []
        for line in inf:
            ll = line.split()
            barcode_seq = ll[0]

            BCnum += 1
            if BCnum <= topN:

                cell_barcode = batch_name+":"+barcode_seq
                if cell_barcode_dict.has_key(cell_barcode):
                    newll = ll + [cell_barcode_dict[cell_barcode]]
                    real_barcodes.append(barcode_seq)
                    cell_vs_adjbarcode_dict[cell_barcode_dict[cell_barcode]] = [batch_name, barcode_seq]
                else:
                    newll = ll + ["NA"]
#                newll += ["NA"]*topN
                #print batch_name,newll
        #print batch_name,real_barcodes
            else:
                p1 = single_barcode_adj(barcode_seq.split("_")[0],5)
                p2 = single_barcode_adj(barcode_seq.split("_")[1],7)
                adj_cell_barcode = batch_name + ":" + p1 + "_" + p2
                if cell_barcode_dict.has_key(adj_cell_barcode):
                    adj_cellname = cell_barcode_dict[adj_cell_barcode]+"_ex"
                    newll = ll + [adj_cellname]
                    if not cell_vs_adjbarcode_dict.has_key(adj_cellname):
                        cell_vs_adjbarcode_dict[adj_cellname] = [batch_name, barcode_seq]
                    else:
                        cell_vs_adjbarcode_dict[adj_cellname].append(barcode_seq)
                else:
                    newll = ll + ["NA"]
 #               for i in range(topN):
 #                   dis_real_barcode.append(ED_twopair(barcode_seq,real_barcodes[i]))
            #newll = ll + ["NA"] + dis_real_barcode
            outf.write("\t".join(map(str,newll))+"\n")
        outf.close()

cell_barcode_outf = open("cell_barcode_list.txt",'w')
for cellname in sorted(cell_vs_adjbarcode_dict.keys()):
    newll = [cellname] + cell_vs_adjbarcode_dict[cellname]
    cell_barcode_outf.write("\t".join(newll)+"\n")
cell_barcode_outf.close()






