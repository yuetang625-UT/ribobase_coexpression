######version 3 week3 spearman correlation for certain_id list
######raw, dedup, all correlation
######gene filter
######run under python manage.py runscript with parameters
######this script is only suit for human data

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

####parameter check
def current_parameters(ribinput, GSMinput, mmin, mmax, outdir):
    #print(ribinput, GSMinput, mmin, mmax, outdir)
    print("input study data folder: ",ribinput)
    print("input samples list: ", GSMinput)
    print("lower ribo length for summary: ", mmin)
    print("upper ribo length for summary: ", mmax)
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
def combine_ribo(GSMinput,ribinput,mmin,mmax,outdir,status):
    file_dict=study_folder(GSMinput)
    new=pd.DataFrame()
    for k,v in file_dict.items():
        try:
            os.chdir(ribinput+"%s"%(k)+"/ribo/experiments")
        except FileNotFoundError as e:
            print(e)
        for j in v:
            try:
                ribo_data("%s.ribo"%j)
            except FileNotFoundError as e:
                print(e)
            else:
                try:
                    CDS_count_length_sum=CDS_count(mmin,mmax,j,status)
                except IndexError as e:
                    print(e)
                except ribopy.core.exceptions.AliasError as e:
                    print(e)
                else:
                    try:
                        new=pd.merge(new,CDS_count_length_sum,on='transcript')
                    except KeyError:
                        new=pd.DataFrame(index=CDS_count_length_sum.index)                                        
    display_correlation(new,outdir,status)
    return new

def run(*args):
    ####parameter setting
    """
    usage: python manage.py runscript my script --script-args="--mmin=XXX"
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--ribinput", help='all analysis folder',default='/scratch/users/yliu5/sample_ribo_data/')
    parser.add_argument("--GSMinput", help='input experiment files in csv', default='/scratch/users/yliu5/ribo_cancer/input.csv')
    parser.add_argument("--mmin", help='lower ribo length', default=28)
    parser.add_argument("--mmax", help='upper ribo length', default=36)
    parser.add_argument("--outdir", help='output directory', default='/scratch/users/yliu5/ribo_cor_analysis/3_15_22/')
    try:
        split_args = shlex.split(args[0])
    except IndexError:
        split_args = []
    args2 = parser.parse_args(split_args)
    current_parameters(args2.ribinput, args2.GSMinput, args2.mmin, args2.mmax, args2.outdir)
    
    #####run main script
    print("#################")
    print("raw")
    raw=combine_ribo(args2.GSMinput,args2.ribinput,args2.mmin,args2.mmax,args2.outdir,status="raw")
    print("dedup")
    dedup=combine_ribo(args2.GSMinput,args2.ribinput,args2.mmin,args2.mmax,args2.outdir,status="dedup")
    inter_overall=pd.merge(raw,dedup,on='transcript')
    display_correlation(inter_overall,args2.outdir,"all") 
