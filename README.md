![QBRC_logo](https://github.com/jcao89757/SCINA/blob/master/QBRC.jpg)
# Tessa
## Introduction
Tessa is a Bayesian model to intergrate TCR sequence profiling with transcriptomes of T cells. Enabled by the recently developed single cell sequencing techniques, which provide both TCR sequences and RNA sequences of each T cell concurrently, Tessa maps the functional landscape of the TCR repertoire, and generates insights into understanding human immune response to diseases. As the first part of tessa, BriseisEncoder is employed prior to the Bayesian algorithm to capture the TCR sequence features and create numerical embeddings. Please refer to our paper for more details: ['Mapping the Functional Landscape of TCR Repertoire'](Pending URL), Zhang Z, Xiong D, et al., 2020
##  Instructions
### Dependencies
Python (version 3 or later)

**Python Packges**

R (version 3.5.1 preferred)

**R Packages**
## Guided tutorial 1: 
### Input data
### Suggested pre-processing workflow
### Parameters setting
### TCR network construction
### Example of results
## Guided tutorial 2: 
## Guided tutorial 3:
## Version update
```{Shell}
python3 BriseisEncoder.py -tcr TestCDR3.csv -model Trained_encoder.h5 -embeding_vectors Atchley_factors.csv -output_TCR test.csv -output_VJ testVJ.csv -output_log test.log
```
