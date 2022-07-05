nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_D_sc_19 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_D_sc_20 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_D_sc_23 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_D_sc_28 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_D_sc_29 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_Dpool_12 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_Dpool_15 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_Dpool_16 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT6_D_bulk_25 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT6_D_bulk_27 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_Dpool_18 &
nohup python ~/Main/Script/single_cell/Tanglab/summaryBarcode_ATACfq.py -i CAT3_D_sc_21 &

wait;

python match_meta_ED.py

python ../ATAC/split_ATACfq_withBCex.py

