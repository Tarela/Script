import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"

#gained_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/fz_combine_all_features/f1_combined_csv/T-ALL_CTCF_gained_AllFeatures.csv'
#lost_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/fz_combine_all_features/f1_combined_csv/T-ALL_CTCF_lost_AllFeatures.csv'

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


##########################################
## basic models to get basic 
##########################################
#chrom_size_file = '/nv/vol190/zanglab/zw5j/data/Genome/hg38/hg38_clean.chrom.sizes'
#chrom_size_df = pd.read_csv(chrom_size_file,sep='\t',header=None,index_col=0);
#chrom_size_df.columns = ['len']

union_CTCF_kept_ids_sample_thre_2_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/fz_gene_CTCF_mete_info/commonly_used_modules_files/union_CTCF_kept_ids_sample_thre_2.csv'
all_binding_middle_domain_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/fz_gene_CTCF_mete_info/commonly_used_modules_files/all_CTCF_domainInfo.csv'
with open(union_CTCF_kept_ids_sample_thre_2_file) as filter_inf, open(all_binding_middle_domain_file) as all_inf:
    kept_ids = [int(i.strip()) for i in filter_inf.readlines()]
    all_bindings = pd.read_csv(all_inf,sep='\t')
all_bindings = all_bindings.loc[kept_ids]
#print(all_bindings)
#exit()
#
gained_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/fz_combine_all_features/f1_combined_csv/T-ALL_CTCF_gained_AllFeatures.csv'
lost_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/fz_combine_all_features/f1_combined_csv/T-ALL_CTCF_lost_AllFeatures.csv'
constitutive_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f4_lost_on_consistent_binding/f1_consistent_bindings/f1_consistent_bindings/constitutive_CTCF_bindings_thre0.8.bed'
#
with open(gained_file) as gained_inf, open(lost_file) as lost_inf, open(constitutive_file) as constitutive_inf:
    gained_df = pd.read_csv(gained_inf,index_col=3)
    lost_df = pd.read_csv(lost_inf,index_col=3)
    const_df = pd.read_csv(constitutive_inf,sep='\t',index_col = 3,header=None)
    const_df.columns = ['chr','start','end','score','strand']

##########################################
##########################################


def mirror_concat(df):

    for column_pos in np.arange(1,len(df.columns)):
        column = df.columns[column_pos]
        mir_column = '-{}'.format(column)
        new_values = np.append([0]*column_pos,df[column].values)[:len(df.index)]
        df[mir_column]=new_values
        #print(df[column],new_values);exit()
    return df
    

def write_out_binding_interactions_all_chroms(binding_df,data,normalization,viewregion,resolution,flag,outdir):
    
    binding_data_out = open(outdir+os.sep+'{}_{}_res{}_{}.csv'.format(data,normalization,resolution,flag),'w')
    
    columns = np.arange(-1*viewregion+resolution,viewregion,resolution)
    binding_data_out.write('{}\t{}\n'.format('id','\t'.join(map(str,columns))))
            
    for chrom in chroms:
        #chrom='chr22'
        print(chrom)
        matrix_file = 'f1_viewpoint_2M_interaction_abstraction/{}_{}_{}_res_{}_view_region_2000000.csv'.format(data,normalization,chrom,resolution)
        if os.path.isfile(matrix_file):
            with open(matrix_file) as matrix_inf:
                matrix_df = pd.read_csv(matrix_inf,sep='\t',index_col=0)
        #####
        # process the matrix df
        #####
        
            # normalization of the data, each interaction score is divided by average interaction score with same distances
            matrix_df = matrix_df/matrix_df.mean()#;print(matrix_file)
            matrix_df = mirror_concat(matrix_df)#;print(binding_data_df.loc[248900000][map(str,columns)].values);exit()
            binding_df_chr = binding_df.loc[binding_df['chr']==chrom]
            for binding_id in binding_df_chr.index:
                binding_pos = all_bindings.loc[binding_id].middle
                view_pos = binding_pos//resolution*resolution
                view_pos_binding_data = matrix_df.loc[view_pos][map(str,columns)].values  
                view_pos_binding_data = np.round(view_pos_binding_data,2)
                #print('\t'.join(map(str,view_pos_binding_data)))
                #print(view_pos_binding_data,len(columns),len(view_pos_binding_data));exit()
                binding_data_out.write('{}\t{}\n'.format(binding_id,'\t'.join(map(str,view_pos_binding_data))))
                #print(view_pos_binding_data);exit()
    binding_data_out.close() #;exit() 
     
     
     
    
        
        
def main(view_region,normalization,chrom):
    
    outdir = 'f2_union_binding_view2M_bothside_interaction'
    os.makedirs(outdir,exist_ok=True)
    
    viewregion = 2000000
    
    for normalization in ['raw','iced']:
        for data in ['A549','HCT116','lung','trans_colon1','trans_colon2']:
            for resolution in [5000,10000,25000]:          
                write_out_binding_interactions_all_chroms(all_bindings,data,normalization,viewregion,resolution,'union',outdir)







if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--viewregion', action = 'store', type = int,dest = 'viewregion', help = 'input file of', metavar = '<int>')
    parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.viewregion,args.normalization,args.chrom)
