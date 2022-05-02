######version 8 spearman correlation for certain_id list
######raw, dedup, all correlation
######CDS percentage
######run under python manage.py runscript with parameters
######this script is only suit for human data
######add log2(upper-quartile normalized + 1)
######fix bug for raw/dedup
######CDS reads vs. length distribution
######bug: first sample not in final file

##############dynamic cutoff for min & max length

import ribopy
from ribopy import Ribo
import matplotlib.pyplot as plt
import seaborn as sns
import os
from metadata.models import *
import pandas as pd
import django
import argparse
import shlex
import numpy as np
import math
import scipy.stats as st
from bioinfokit.analys import norm


####parameter check
def current_parameters(ribinput, GSMinput, outdir):
    #print(ribinput, GSMinput, mmin, mmax, outdir)
    print("input study data folder: ",ribinput)
    print("input samples list: ", GSMinput)
    print("output results: ", outdir)
    
###read input experiment
def read_table(Experiment_file):
    table=pd.read_csv(Experiment_file)
    return table

####study_id and and experiment link
def study_folder(Experiment_file):
    study_dic={}
    experiment_list=read_table(Experiment_file)['experiment_alias'].values.tolist()
    for i in experiment_list:
        study=Experiment.objects.get(experiment_alias=i).study
        if study not in study_dic:
            study_dic[study]=[]
            study_dic[study].append(i)
        else:
            study_dic[study].append(i)
    return study_dic

####read human gene name
def defalternative_human_alias(x):
    x_pieces = x.split("|")
    return x_pieces[1] + "_" + x_pieces[5]

####read ribo files
def ribo_data(ribo_path):
    ribo_object = Ribo(ribo_path,alias = defalternative_human_alias)
    return ribo_object

####read count file by file    
def CDS_count(mmin,mmax,experiment_id,status):
    ribo_object = ribo_data("%s.ribo"%experiment_id)
    #ribo_object.print_info() #print the information
    CDS_count_length_sum = ribo_object.get_region_counts(region_name = "CDS",
                              range_lower    = int(mmin),
                              range_upper    = int(mmax),
                              sum_lengths    = True,
                              sum_references = False,
                              alias          = True) ### get data sum across
    CDS_count_length_sum_fin=CDS_count_length_sum.add_suffix(status)
    return CDS_count_length_sum_fin

####reads CDS percentage all not length 
'''def CDS_percentage(experiment_id,mmax,mmin):
    ribo_object = ribo_data("%s.ribo"%experiment_id)
    total_count = ribo_object.info['experiments']['%s'%experiment_id]['Reads']
    CDS_counts = ribo_object.get_region_counts(region_name='CDS')
    CDS_per_fin = CDS_counts/total_count
    return CDS_per_fin'''

#####dynamic cutoff
def intevl(experiment_id):
    ribo_object = ribo_data("%s.ribo"%experiment_id)
    data=ribo_object.get_length_dist("CDS")
    data.reset_index(inplace=True)
    pct_85=sum(data["%s"%experiment_id])*0.85
    pct_90=sum(data["%s"%experiment_id])*0.90
    value=data[data["%s"%experiment_id]==data["%s"%experiment_id].max()]["%s"%experiment_id].values[0]
    mmin=mmax=data[data["%s"%experiment_id]==data["%s"%experiment_id].max()]['read_length'].values[0]
    if mmax<23:
        mmin=mmax="warming:short peak"
    else:
        while value<=pct_85 :
            if mmax <40 and mmin>20 :
                if data[data['read_length']==mmax+1]["%s"%experiment_id].values[0] >= data[data['read_length']==mmin-1]["%s"%experiment_id].values[0]:
                    mmax+=1
                    value+=data[data['read_length']==mmax]["%s"%experiment_id].values[0]
                else:
                    mmin-=1
                    value+=data[data['read_length']==mmin]["%s"%experiment_id].values[0]
            elif mmax==40 and mmin>20:
                mmin-=1
                value+=data[data['read_length']==mmin]["%s"%experiment_id].values[0]
            elif mmax <40 and mmin==20:
                mmax+=1
                value+=data[data['read_length']==mmax]["%s"%experiment_id].values[0]
            else:
                break
    #print(min,max)
    read_pct=value/sum(data["%s"%experiment_id])
    return mmin,mmax,read_pct
    
####CDS coverage and percentage
def CDS_coverage_and_pct(experiment_id,mmin,mmax):
    ribo_object = ribo_data("%s.ribo"%experiment_id)
    df=ribo_object.get_length_dist("CDS")
    df.reset_index(inplace=True)
    df_keep=df[(df['read_length']>(mmin-1)) & (df['read_length']<(mmax+1))]
    cov=sum(df_keep['read_length']*df_keep["%s"%experiment_id])/33636272
    #read_pct=sum(df_keep['read_length']*df_keep["%s"%experiment_id])/sum(df['read_length']*df["%s"%experiment_id])
    #summ=df.values
    CDS_total_count=df_keep['%s'%experiment_id].sum()
    UTR3_counts = ribo_object.get_region_counts(region_name='UTR3',range_lower = int(mmin),range_upper = int(mmax))
    UTR5_counts = ribo_object.get_region_counts(region_name='UTR5',range_lower = int(mmin),range_upper = int(mmax))
    UTR3_j_counts = ribo_object.get_region_counts(region_name='UTR3_junction',range_lower = int(mmin),range_upper = int(mmax))
    UTR5_j_counts = ribo_object.get_region_counts(region_name='UTR5_junction',range_lower = int(mmin),range_upper = int(mmax))
    CDS_pct=CDS_total_count/(CDS_total_count+UTR3_counts+UTR5_counts+UTR3_j_counts+UTR5_j_counts)
    return cov,CDS_pct

#####quartile normalized
def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    df_log=np.log2(df_qn+1)
    return(df_log)

def CPM_normalize(df):
    # now, normalize raw counts using CPM method 
    nm = norm()
    nm.cpm(df=df)
    # get CPM normalized dataframe
    cpm_df = nm.cpm_norm
    df_log=np.log2(cpm_df+1)
    return df_log

######gene qualitiy controls
def drop_almost_zero(df, percentage):
    row_cut_off = int(percentage/100*len(df.columns))
    df = df[(df==0).sum(axis='columns') <= row_cut_off]

    column_cut_off = int(percentage/100*len(df)) 
    b = (df == 0).sum(axis='rows')
    df = df[ b[ b <= column_cut_off].index.values ]

    return df

######gene level correlation
def gene_cor(df_log):
    df_log_t=df_log.T
    return df_log_t

####correlation function
def display_correlation(df,path,status):
    r = df.corr(method="spearman")
    fig=plt.figure(figsize=(11,11))
    heatmap = sns.heatmap(df.corr(), vmin=0, 
                      vmax=1, annot=True)
    plt.title("Spearman Correlation")
    fig.savefig(path+'%s_cor.png'%(status), dpi=fig.dpi)
    plt.close("all")
    r.to_csv(path+'%s_cor.csv'%(status),index=True)
    return r

####combine individual dataframe and caculate cor
def combine_ribo(GSMinput,ribinput,outdir,status):
    CDS_per={}
    CDS_region={}
    file_dict=study_folder(GSMinput)
    new=pd.DataFrame()
    for k,v in file_dict.items():
        if status == "dedup":
            data_dir=ribinput+"%s_dedup"%(k)+"/ribo/experiments"
        else:
            data_dir=ribinput+"%s"%(k)+"/ribo/experiments"
        try:
            os.chdir(data_dir)
        except FileNotFoundError as e:
            print(e)
            print(k)
        else:
            for j in v:
                try:
                    ribo_data("%s.ribo"%j)
                except FileNotFoundError as e:
                    print(e)
                    print(j)
                else:
                    try:
                        mmin,mmax,select_pct=intevl(j)
                        CDS_count_length_sum=CDS_count(mmin,mmax,j,status)
                    except IndexError as e:
                        print(e)
                    except ribopy.core.exceptions.AliasError as e:
                        print(e)
                    except ValueError as e:
                        print(e)
                    else:
                        try:
                            new=pd.merge(new,CDS_count_length_sum,on='transcript')
                            CDS_per[j]=[CDS_coverage_and_pct(j,mmin,mmax)]
                            CDS_region[j]=[mmin,mmax,select_pct]
                        except KeyError:
                            new=pd.DataFrame(index=CDS_count_length_sum.index)
    #display_correlation(new,outdir,status)
    return new,CDS_per,CDS_region

def run(*args):
    ####parameter setting
    """
    usage: python manage.py runscript my script --script-args="--mmin=XXX"
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--ribinput", help='all analysis folder',default='/scratch/users/mjgeng/process-multiple-ribo/output/')
    parser.add_argument("--GSMinput", help='input experiment files in csv', default='/scratch/users/yliu5/ribo_cancer/input.csv')
    parser.add_argument("--outdir", help='output directory', default='/scratch/users/yliu5/ribo_cor_analysis/4_28_22/')
    try:
        split_args = shlex.split(args[0])
    except IndexError:
        split_args = []
    args2 = parser.parse_args(split_args)
    current_parameters(args2.ribinput, args2.GSMinput, args2.outdir)
    
    #####run main script
    print("#################")
    #print("raw")
    #raw,raw_CDS_percentage=combine_ribo(args2.GSMinput,args2.ribinput,args2.mmin,args2.mmax,args2.outdir,status="raw")
    #raw_per_df=pd.DataFrame(raw_CDS_percentage.items())
    #raw.to_csv("/scratch/users/yliu5/ribo_cor_analysis/4_18_22/cancer_ribo_raw.csv")
    #raw_per_df.to_csv("/scratch/users/yliu5/ribo_cor_analysis/4_20_22/cancer_cdscov_raw.csv")
    print("dedup")
    ribo_dedup,dedup_CDS_percentage,dedup_CDS_region=combine_ribo(args2.GSMinput,args2.ribinput,args2.outdir,status="dedup")
    dedup_per_df=pd.DataFrame(dedup_CDS_percentage.items())
    dedup_region=pd.DataFrame(dedup_CDS_region.items())
    ribo_dedup.to_csv("/scratch/users/yliu5/ribo_cor_analysis/4_28_22/cancer_ribo_dedup.csv")
    dedup_per_df.to_csv("/scratch/users/yliu5/ribo_cor_analysis/4_28_22/cancer_readpct_dedup.csv")
    dedup_region.to_csv("/scratch/users/yliu5/ribo_cor_analysis/4_28_22/cancer_readregion_dedup.csv")
    
    

