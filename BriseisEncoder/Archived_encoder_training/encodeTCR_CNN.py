import numpy as np
import pandas as pd
import os
from keras.layers import Input, Conv2D,AveragePooling2D, UpSampling2D, Dropout, Flatten, Dense, Reshape,BatchNormalization
from keras.models import Model, Sequential, load_model
from keras.optimizers import Adam
import sys
sys.path.append('/project/bioinformatics/Xiao_lab/s421955/projects/scTCR/code')
import csv
import matplotlib.pyplot as plt
from encodeTCR_functions import preprocessMultiInputs, preprocess, datasetMap,aamapping,plotLoss,plotHeat
import random as rd
#args = sys.argv
#tcr_dir = args[args.index('-tcr') + 1]
#1. read TCR contigs
prediction_tcr_dir='/home2/s421955/projects/scTCR/data/GEX_VDJ_SampleT/outs/filtered_contig_annotations.csv'
tcr_dir='/project/bioinformatics/Xiao_lab/s421955/projects/scTCR/data/autoencoder_index/TCRfilesConfig_pretrain_random.txt'
#test=preprocessMultiInputs(tcr_dir,'test')
train=preprocessMultiInputs(tcr_dir,'train')
k=round(train.shape[0]*0.1)
test_ind=rd.sample(range(0,train.shape[0]),k=k)
test=train.loc[test_ind]
train=train.loc[~train.index.isin(test_ind)]
train.index=range(0,train.shape[0])
test.index=range(0,test.shape[0])
tcr=preprocess(prediction_tcr_dir)
#2. read Factor analysis transition mat for TCR contigs, analysis from
# [https://www.ncbi.nlm.nih.gov/pubmed/15851683]
aa_idx_dir='/project/bioinformatics/Xiao_lab/s421955/projects/scTCR/data/factor_analysis_index.csv'
aa_dict=dict()
with open(aa_idx_dir,'r') as aa:
    aa_reader=csv.reader(aa)
    next(aa_reader, None)
    for rows in aa_reader:
        aa_name=rows[0]
        aa_factor=rows[1:len(rows)]
        aa_dict[aa_name]=np.asarray(aa_factor,dtype='float')

#3. map TCR contigs
TCR_dict_train = datasetMap(train,aa_dict)
TCR_dict_test = datasetMap(test,aa_dict)
TCR_dict=datasetMap(tcr,aa_dict)
encoding_dim=80
#6. convert to 3D array
TCR_train=np.stack(TCR_dict_train.values())
TCR_train=TCR_train.reshape(-1,encoding_dim,5,1)
TCR_test=np.stack(TCR_dict_test.values())
TCR_test=TCR_test.reshape(-1,encoding_dim,5,1)
TCR_contigs=np.stack(TCR_dict.values())
TCR_contigs=TCR_contigs.reshape(-1,encoding_dim,5,1)
#5. Build encoders
max_features=30
batch_size=256
input_TCR=Input(shape=(encoding_dim,5,1))
def TCRencoder(input_TCR):
    #encoder
    conv1=Conv2D(30,(5,2),activation='selu')(input_TCR)
    conv1_bn=BatchNormalization(axis=1)(conv1)
    pool1 = AveragePooling2D(pool_size=(4, 1))(conv1_bn)
 
    conv2=Conv2D(20,(4,2),activation='selu')(pool1)
    conv2_bn=BatchNormalization(axis=1)(conv2)
    pool2 = AveragePooling2D(pool_size=(4, 1))(conv2_bn)
 
    pool2_flat=Flatten()(pool2)
    dense1=Dense(30,activation="selu")(pool2_flat)
    dense1_do=Dropout(0.01)(dense1)
 
    #bottle_neck
    bottle_neck=Dense(30,activation="selu")(dense1_do)
    bottle_neck_do=Dropout(0.01)(bottle_neck)
 
    #decoder
    dense2=Dense(30,activation="selu")(bottle_neck_do)
    dense2_do=Dropout(0.01)(dense2)
 
    dense3=Dense(240,activation="selu")(dense2_do)
    dense3_3d=Reshape((4,3,20))(dense3)
    dense3_3d_bn=BatchNormalization(axis=1)(dense3_3d)
 
    up1=UpSampling2D((5,2))(dense3_3d_bn)
    conv3=Conv2D(30,(4,3),activation='selu')(up1)
    conv3_bn=BatchNormalization(axis=1)(conv3)
 
    up2=UpSampling2D((5,2))(conv3_bn)
    conv4=Conv2D(1,(6,4),activation='linear')(up2)
 
    return conv4
 
TCRencoder=Model(input_TCR,TCRencoder(input_TCR))
TCRencoder.compile(optimizer=Adam(), loss='mean_squared_error')
print(TCRencoder.summary())
history=TCRencoder.fit(TCR_train, TCR_train,epochs=50,batch_size=batch_size,shuffle=True,validation_data=(TCR_test, TCR_test))
#TCRencoder.save('/project/bioinformatics/Xiao_lab/s421955/projects/scTCR/data/TCRencodermodel_121118.h5')
TCRencoder=load_model('/project/bioinformatics/Xiao_lab/s421955/projects/scTCR/data/TCRencodermodel_121118.h5')
encoder=Model(TCRencoder.input,TCRencoder.layers[-12].output)
encoded_mat=encoder.predict(TCR_contigs)
encoded_mat=pd.DataFrame(encoded_mat,index=tcr['cdr3'])
encoded_mat.to_csv('/home2/s421955/projects/scTCR/data/GEX_VDJ_SampleT/outs/encoded_filtered_contig_pretrain.csv',sep=',')
#plot loss
# summarize history for loss
plotLoss(history)
#draw heatmaps to test decoder accuracy
decoded_mat=TCRencoder.predict(TCR_contigs)
plotHeat(TCR_contigs[0],decoded_mat[0],'Original','Decoded')
plotHeat(TCR_contigs[1],decoded_mat[1],'Original','Decoded')
plotHeat(TCR_contigs[2],decoded_mat[2],'Original','Decoded')
plotHeat(TCR_contigs[3],decoded_mat[3],'Original','Decoded')
plotHeat(TCR_contigs[4],decoded_mat[4],'Original','Decoded')
