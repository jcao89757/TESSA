#Parameters:
#tcr_dir: a .csv file to input the CDR3 sequences, contains at least 2 columns naming 'contig_id' and 'cdr3'. Two optional
# 		  columns 'v_gene' and 'j_gene' denotes the V gene and J gene subgroups (TRBV1-30 and TRBJ1/2) of TRB recombinants.
#		  'contig_id' should contains non-repeated ID strings, and 'cdr3' should contains valid TRB CDR3 peptide sequences.
#model_dir: a .h5 file containing trained encoding models from CDR3 seqs called from bulk-tumor RNA-seq data.
#aa_dict_dir: a .csv file storing Atchley's amino acid embedding, details in https://www.ncbi.nlm.nih.gov/pubmed/15851683.
#output_embeddedTCR_dir: a .json file dir to store Atchley factor encoded CDR3 peptide seqs.
#output_decodedTCR_dir: a .json file dir to store deencoded CDR3 peptide seqs.
#output_log_dir: a plain text log-file to record any errors or warnings from this script.
##########################################################################################################################
#Example:
#python3 OutputAccTest.py -tcr ../example_data/example_TCRmeta.csv -model TrainedEncoder.h5 -embeding_vectors Atchley_factors.csv \
#-embeddedTCR_dir test1.json -decodedTCR_dir test2.json -output_log test.log
##########################################################################################################################
#Import dependencies
import sys
sys.path.append('/home2/s421955/.conda/envs/keras_test/site-packages')
import numpy as np
import pandas as pd
import os
import csv
from keras.models import load_model
from keras.models import Model
import json
#Read data
args = sys.argv
tcr_dir = args[args.index('-tcr')+1]
model_dir=args[args.index('-model')+1]
aa_dict_dir=args[args.index('-embeding_vectors')+1]
output_embeddeedTCR_dir=args[args.index('-embeddedTCR_dir')+1]
output_decodedTCR_dir=args[args.index('-decodedTCR_dir')+1]
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
    })
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
encoder=Model(TCRencoder.input,TCRencoder.output)
encoded_mat=encoder.predict(TCR_contigs)
with open(output_embeddeedTCR_dir,'w+') as f:
    json.dump(TCR_contigs.tolist(),f)
f.close()
with open(output_decodedTCR_dir,'w+') as f:
    json.dump(encoded_mat.tolist(),f)
f.close()
print('Mission Accomplished.\n')
log_file.close()
