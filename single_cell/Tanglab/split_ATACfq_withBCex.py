import os,sys,gzip
summary_folder = "/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/fastq/05_20181030_mESC_tagmentation_CAbeads/split_barcode_fastq/summary/"
data_folder = "/nm/vol190/zanglab/sh8tv/Project/scATAC/Data/scATAC/Tanglab_scATAC/fastq/05_20181030_mESC_tagmentation_CAbeads/combined_ATAC/"

this_cell = sys.argv[1]

startpos = 27
barcodeLen = 8
inf = open(summary_folder+"cell_barcode_list.txt")
for line in inf:
    ll = line.split()
    cell_name = ll[0] 
    batch_name = ll[1]
    usebarcodes = ll[2:]

    if not cell_name == this_cell:
        continue

    fq1 = "%s/%s/%s_R1.fastq.gz"%(data_folder,batch_name,batch_name)
    fq2 = "%s/%s/%s_R2.fastq.gz"%(data_folder,batch_name,batch_name)
    inf1 = gzip.open(fq1,"rb")
    inf2 = gzip.open(fq2,"rb")

    line_count = 0

    outf1 = open("%s_R1.fastq"%(cell_name),"w")
    outf2 = open("%s_R2.fastq"%(cell_name),"w")

    while 1:
        line_count += 1
        p1 = inf1.readline()
        p2 = inf2.readline()

        if line_count%4 == 1:
            p1l1 = p1
            p2l1 = p2
        elif line_count%4 == 2:
            p1l2 = p1[startpos:]
            p2l2 = p2[startpos:]
            line_p5 = p1[:barcodeLen]
            line_p7 = p2[:barcodeLen] 
            #print "p1",p1,line_p5
            #print "p2",p2,line_p7               
        elif line_count%4 == 3:
            p1l3 = p1
            p2l3 = p2
        else:
            p1l4 = p1[startpos:]
            p2l4 = p2[startpos:]

            barcode = line_p5 + "_" + line_p7
            if barcode in usebarcodes:
                outf1.write(p1l1)
                outf1.write(p1l2)
                outf1.write(p1l3)
                outf1.write(p1l4)
                outf2.write(p2l1)
                outf2.write(p2l2)
                outf2.write(p2l3)
                outf2.write(p2l4)

        if p1.strip() == "" and p2.strip() == "":
            break
    inf1.close()
    inf2.close()
    outf1.close()
    outf2.close()








