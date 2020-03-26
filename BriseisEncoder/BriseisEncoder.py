##########################################################################################################################
#BriseisEncoder 1.0.0: A deep learning encoder for CDR3 sequences in the TCR space
##########################################################################################################################
#BriseisEncoder is capable of transforming CDR3 peptide sequences from productive
#TCR-beta chains into 15-digit informative numerical vectors.
#01312019 Ze Zhang <Ze.Zhang@UTsouthwestern.edu>
##########################################################################################################################
#Dependencies: 
#Python 3.6.4 (preferred)
#numpy (>=1.15.4), pandas (>=0.23.4), keras (>=2.2.4)
##########################################################################################################################
#Parameters:
#tcr_dir: a .csv file to input the CDR3 sequences, contains at least 2 columns naming 'contig_id' and 'cdr3'. Two optional
# 		  columns 'v_gene' and 'j_gene' denotes the V gene and J gene subgroups (TRBV1-30 and TRBJ1/2) of TRB recombinants.
#		  'contig_id' should contains non-repeated ID strings, and 'cdr3' should contains valid TRB CDR3 peptide sequences.
#model_dir: a .h5 file containing trained encoding models from CDR3 seqs called from bulk-tumor RNA-seq data.
#aa_dict_dir: a .csv file storing Atchley's amino acid embedding, details in https://www.ncbi.nlm.nih.gov/pubmed/15851683.
#output_encodedTCR_dir: a .csv file dir to store encoded CDR3 peptide seqs.
#output_log_dir: a plain text log-file to record any errors or warnings from this script.
#output_encodedVJ_dir: a .csv file dir to store one-hot encoded TRBV/TRBJ genes.
##########################################################################################################################
#Example:
#python3 BriseisEncoder.py -tcr TestCDR3.csv -model Trained_encoder.h5 -embeding_vectors Atchley_factors.csv \
#-output_TCR test.csv -output_VJ testVJ.csv -output_log test.log
##########################################################################################################################
#Import dependencies
import numpy as np
import pandas as pd
import os
import csv
import sys
from keras.models import load_model
from keras import Model
#Read data
args = sys.argv
tcr_dir = args[args.index('-tcr')+1]
model_dir=args[args.index('-model')+1]
aa_dict_dir=args[args.index('-embeding_vectors')+1]
output_encodedTCR_dir=args[args.index('-output_TCR')+1]
if '-output_VJ' in args:
    output_encodedVJ_dir=args[args.index('-output_VJ')+1]
output_log_dir=args[args.index('-output_log')+1]
encode_dim=80
#Define functions
def preprocess(filedir):
    #Preprocess TCR files
    print('Processing: '+filedir)
    if not os.path.exists(filedir):
        print('Invalid file path: ' + filedir)
        return 0
    dataset = pd.read_csv(filedir, header=0)
    if dataset.isnull().values.any():
        print('Input data contains NAs.')
        dataset=dataset.dropna()
    data_new=pd.DataFrame({
        'contig_id':dataset['contig_id'],
        'cdr3':dataset['cdr3'],
        'v_gene':None,
        'j_gene':None
    })
    if 'v_gene' in dataset.columns:
        data_new['v_gene']= dataset['v_gene']
    else:
        data_new=data_new.drop(columns='v_gene')
    if 'j_gene' in dataset.columns:
        data_new['j_gene'] = dataset['j_gene']
    else:
        data_new=data_new.drop(columns='j_gene')
    data_new.index=range(0,dataset.shape[0])
    return data_new

def aamapping(peptideSeq,aa_dict,encode_dim):
    #Transform aa seqs to Atchley's factors.
    peptideArray = []
    if len(peptideSeq)>encode_dim:
        print('Length: '+str(len(peptideSeq))+' over bound!')
        peptideSeq=peptideSeq[0:encode_dim]
    for aa_single in peptideSeq:
        try:
            peptideArray.append(aa_dict[aa_single])
        except KeyError:
            print('Not proper aaSeqs: '+peptideSeq)
            peptideArray.append(np.zeros(5,dtype='float64'))
    for i in range(0,encode_dim-len(peptideSeq)):
        peptideArray.append(np.zeros(5,dtype='float64'))
    return np.asarray(peptideArray)

def datasetMap(dataset,aa_dict,encode_dim):
    #Wrapper of aamapping
    TCR_dict=dict()
    for i in range(0,len(dataset['cdr3'])):
        TCR_key=dataset['contig_id'][i]
        if TCR_key in TCR_dict.keys():
            TCR_key=TCR_key+'_'+str(i)
        TCR_dictarray=aamapping(dataset['cdr3'][i],aa_dict,encode_dim)
        TCR_dict[TCR_key]=TCR_dictarray
    return TCR_dict

def embedVJ(genelist,maplist):
    #Embed VJgenes
    VJ_array=[]
    for gene in genelist:
        ind=np.zeros(len(maplist))
        try:
            find=maplist.index(gene)
            ind[find]=1
            VJ_array.append(ind)
        except ValueError:
            print('Gene out of bound!'+gene)
            VJ_array.append(ind)
            next
    return np.asarray(VJ_array)

#Main functions, data preprocess
log_file=open(output_log_dir,'w')
sys.stdout=log_file
print('Mission loading.')
tcr=preprocess(tcr_dir)
aa_dict=dict()
with open(aa_dict_dir,'r') as aa:
    aa_reader=csv.reader(aa)
    next(aa_reader, None)
    for rows in aa_reader:
        aa_name=rows[0]
        aa_factor=rows[1:len(rows)]
        aa_dict[aa_name]=np.asarray(aa_factor,dtype='float')
TCR_dict=datasetMap(tcr,aa_dict,encode_dim)
TCR_contigs=np.stack(TCR_dict.values())
TCR_contigs=TCR_contigs.reshape(-1,encode_dim,5,1)
#Model prediction
TCRencoder=load_model(model_dir)
encoder=Model(TCRencoder.input,TCRencoder.layers[-12].output)
encoded_mat=encoder.predict(TCR_contigs)
encoded_mat=pd.DataFrame(encoded_mat,index=tcr['contig_id'])
encoded_mat.to_csv(output_encodedTCR_dir,sep=',')
if 'v_gene' in tcr.columns or 'j_gene' in tcr.columns:
    maplist=['TRBV'+str(i) for i in range(1,31)]
    maplist.append('TRBJ1')
    maplist.append('TRBJ2')
    if 'v_gene' in tcr.columns:
        v_map=embedVJ(tcr['v_gene'],maplist)[:,0:30]
    else:
        print('V genes are missing!')
        v_map=np.zeros((tcr.shape[0],30),dtype='float64')
    if 'j_gene' in tcr.columns:
        j_map = embedVJ(tcr['j_gene'], maplist)[:, 30:32]
    else:
        print('J genes are missing!')
        j_map=np.zeros((tcr.shape[0],2),dtype='float64')
    VJ_map=np.concatenate((v_map,j_map),axis=1)
    VJ_map=pd.DataFrame(data=VJ_map,index=tcr['contig_id'],columns=maplist)
    VJ_map.to_csv(output_encodedVJ_dir,sep=',')
print('Mission Accomplished.\n')
log_file.close()
