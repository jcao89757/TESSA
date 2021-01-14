#tcr_dir inputs a txt file, each line contains either a training table (with path, csv) or a testing table (with path)
#files should all have the same format with cellRanger outputs. minimum 4 columns ['contig_id','barcode',chain,cdr3]
#training files are listed between '#train' line and '#trainEOF', testing files under '#test' line and '#testEOF'
import numpy as np
import pandas as pd
import os
import sys
import csv
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.preprocessing import OneHotEncoder

def preprocess(filedir,vjgene=False):
    print(filedir)
    dataset = pd.read_csv(filedir, header=0)
    dataset['contig_id'] = [b.replace('-', '.') for b in dataset['contig_id']]
    dataset['barcode'] = [b.replace('-', '.') for b in dataset['barcode']]
    # filter and keep cell owned, high_confidence TCR contigs
    dataset=dataset.loc[dataset['is_cell']]
    dataset = dataset.loc[dataset['productive']]
    dataset = dataset.loc[dataset['high_confidence']]
    if vjgene:
        if 'label' in dataset.columns:
            dataset = dataset[['barcode', 'contig_id', 'chain', 'label','cdr3', 'v_gene', 'j_gene']]
            dataset = dataset.loc[dataset['label'] != 999]
        else:
            dataset = dataset[['barcode', 'contig_id', 'chain', 'cdr3','v_gene','j_gene']]
    else:
        if 'label' in dataset.columns:
            dataset = dataset[['barcode', 'contig_id', 'chain', 'label','cdr3']]
            dataset = dataset.loc[dataset['label'] != 999]
        else:
            dataset = dataset[['barcode', 'contig_id', 'chain', 'cdr3']]
    dataset = dataset.loc[dataset['cdr3'] != 'None']
    dataset = dataset.loc[dataset['chain']=='TRB']
    dataset.index=range(0,dataset.shape[0])
    return dataset

def preprocessMultiInputs(tcr_dir,train_or_test):
    with open(tcr_dir,'r') as filelist:
        lines=filelist.read().splitlines()
    #find files for 'train' or for 'test'
    filesind=range(lines.index('#'+train_or_test)+1,lines.index('#'+train_or_test+'EOF'))
    print('Number of '+train_or_test+' files: '+str(len(filesind)))
    dataset=pd.DataFrame()
    for i in filesind:
        file=lines[i]
        if not os.path.exists(file):
            print('Invalid file path: '+file)
            break
        onedata=preprocess(file)
        dataset=dataset.append(onedata)
    dataset=dataset.loc[~pd.isnull(dataset['cdr3'])]
    dataset.index = range(0, dataset.shape[0])
    return dataset

def aamapping(peptideSeq,aa_dict,encode_dim):
    #the longest sting will probably be shorter than 80 nt
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

def datasetMap(dataset,aa_dict,encode_dim=30):
    TCR_dict=dict()
    for i in range(0,len(dataset['cdr3'])):
        TCR_key=dataset['contig_id'][i]
        if TCR_key in TCR_dict.keys():
            TCR_key=TCR_key+'_'+str(i)
        TCR_dictarray=aamapping(dataset['cdr3'][i],aa_dict,encode_dim)
        TCR_dict[TCR_key]=TCR_dictarray
    return TCR_dict

def datasetConvert(TCR_dict):
    TCR_array = np.empty((len(TCR_dict),encoding_dim,5),dtype='float64')
    for ik in range(TCR_dict.keys()):
        TCR_array=np.append(TCR_array,TCR_dict[k].reshape(1,encoding_dim*5),axis=0)
        #TCR_array = np.append(TCR_array, TCR_dict[k], axis=0)
    return np.asarray(TCR_array)
def embedVJ(genelist,maplist):
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
#try biovec from https://github.com/kyu999/biovec
#ProtVec: A Continuous Distributed Representation of Biological Sequences

def biovectmap(dataset):
    TCR_array=[]
    for i in range(0,dataset.shape[0]):
        seq=dataset['cdr3'][i]
        TCR_array.append(pv.to_vecs(seq))
    return np.asarray(TCR_array)
def plotLoss(history):
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    if history.history['acc']:
        plt.plot(history.history['acc'])
    if history.history['val_acc']:
        plt.plot(history.history['val_acc'])
    plt.title('model loss')
    plt.ylabel('loss/acc')
    plt.xlabel('epoch')
    plt.legend(['train_loss', 'val_loss','train_acc','val_acc'], loc='upper left')
    plt.show()

def plotHeat(mat2plot1,mat2plot2,title1,title2):
    mat2plot1 = np.reshape(mat2plot1, (80, 5))
    mat2plot1 = mat2plot1.transpose()
    mat2plot2 = np.reshape(mat2plot2, (80, 5))
    mat2plot2 = mat2plot2.transpose()
    fig=plt.figure()
    ax1=fig.add_subplot(2,1,1)
    ax1.set_title(title1)
    ax1.imshow(mat2plot1)
    ax2=fig.add_subplot(2,1,2)
    ax2.set_title(title2)
    ax2.imshow(mat2plot2)
    plt.show()

def plotROC(probs,y_test):
    fpr, tpr, threshold = metrics.roc_curve(y_test, probs)
    roc_auc = metrics.auc(fpr, tpr)
    plt.title('ROC')
    plt.plot(fpr, tpr, 'b', label='AUC = %0.2f' % roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.show()

def checkDictIndex(TCR_dict,dataset):
    TCR_keys=np.stack(TCR_dict.keys())
    if not np.array_equal(TCR_keys[0:1000],dataset['contig_id'][0:1000]):
        print('Unmatched dicts!')
    else:
        print('Check passed!')
