![QBRC_logo](https://github.com/jcao89757/SCINA/blob/master/QBRC.jpg)
# Tessa
## Introduction
Tessa is a Bayesian model to intergrate T cell receptor (TCR) sequence profiling with transcriptomes of T cells. Enabled by the recently developed single cell sequencing techniques, which provide both TCR sequences and RNA sequences of each T cell concurrently, Tessa maps the functional landscape of the TCR repertoire, and generates insights into understanding human immune response to diseases. As the first part of tessa, BriseisEncoder is employed prior to the Bayesian algorithm to capture the TCR sequence features and create numerical embeddings. Please refer to our paper for more details: ['Mapping the Functional Landscape of TCR Repertoire'](Pending URL), Zhang Z, Xiong D, et al., 2020. 

Researchers searching for more bioinformatics tools please visit our lab website: https://qbrc.swmed.edu/labs/wanglab/index.php.
##  Instructions
The tessa algorithm is implemented in python and R. We suggest that users execute the python scripts with Linux shell commands.
### Dependencies
Python (version 3.6.4 preferred), R (version 3.5.1 preferred), Linux (x86_64-redhat-linux-gn) shell (4.2.46(2) preferred)

**Python Packges**

numpy (version 1.15.4 or later), pandas (version 0.23.4 or later), keras (version 2.2.4 or later), os, csv, sys

**R Packages**

Rtsne (version 0.15), MASS (version 7.3-51.4), LaplacesDemon (version 16.1.1) 
## Guided tutorial
In this tutorial, we will show a complete work flow from pre-processing TCR sequences with the BriseisEncoder to constructing TCR networks. The toy example data we used in this tutorial including the TCR sequences and the RNA expression data is availiable [here](https://github.com/jcao89757/TESSA/tree/master/example_data).
### Installation
The full tessa algorithm includes all scripts in the folder [BriseisEncoder](https://github.com/jcao89757/TESSA/tree/master/BriseisEncoder) and [Tessa](https://github.com/jcao89757/TESSA/tree/master/Tessa), and the python script [Tessa_main.py](https://github.com/jcao89757/TESSA/blob/master/Tessa_main.py). Please download all the scripts above, save Tessa_main.py in the same path as the two folders, and keep the directry structure unchanged.
### Input data
The tessa model takes two input data matrices to construct TCR networks. The TCR sequences (or the prepared embedded TCRs, details in **Note 2**), and the single-cell RNA expression levels of the same group of cells.
1. A meta data matrix contains TCR sequences and cell identifiers. TCR sequences used in tessa are the peptide sequences of TCR-beta chain CDR3 regions. Cell identifiers are unique for T cells. They could be self-defined IDs, or cell barcodes, etc. The two columns are required in the meta data matrix, and the column names are specified as 'cdr3' and 'contig_id' **(Fig. 1)**. Each row represent one T cell. Please find the [.csv example ](https://github.com/jcao89757/TESSA/blob/master/example_data/example_TCRmeta.csv) for the meta data matrix. 
2. A matrix representing the gene expression levels. Columns correspond to cells, rows correspond to genes **(Fig. 2)**. Cells should be in the same order as the cell identifiers in the meta data matrix. Please find the the [.csv example](https://github.com/jcao89757/TESSA/blob/master/example_data/example_exp.csv) for the expression matrix.

![meta_example](https://github.com/jcao89757/TESSA/blob/master/example_data/meta_example_fig.png)

**Fig.1 |** An example of a TCR meta data matrix in .csv format.
![exp_example](https://github.com/jcao89757/TESSA/blob/master/example_data/exp_example_fig.png)

**Fig.2 |** An example of an expression data matrix in .csv format.
### Suggested pre-processing workflow
The meta data matrix is suggested to be double-checled and make sure that each element in the column 'contig_id' is unique. The 'cdr3' column allows deplicates, but sequences with any letters, numbers or symbols that do not represent amino acids should be removed. 

For the expression data, we suggest the pre-processing code below to achieve the best performance if the expression matrix compresis raw counts. The log-transformation is always suggested to avoid heavy-tailed datasets.
```{r}
#Install preprocessCore if required
source("http://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
library('preprocessCore')
#Read data
exp_raw=read.csv('your/path/raw_example_exp.csv',row.names=1,stringsAsFactors = F)
#Log scale and quantile normalization
exp_raw=log(exp_raw+1)
exp[]=normalize.quantiles(exp_raw)
```
### Parameters setting

### TCR network construction
### Example of results
## Note 1: Run BriseisEncoder to generate TCR embeddings.

![embedding_example](https://github.com/jcao89757/TESSA/blob/master/example_data/embedding_example_fig.png)

**Fig.3 |** An example of embedded TCRs in .csv format.

## Note 2: Run tessa with prepared embedded TCRs.
## Version update
1.0.0: First release. (03-29-2020)
```{Shell}
python3 BriseisEncoder.py -tcr TestCDR3.csv -model Trained_encoder.h5 -embeding_vectors Atchley_factors.csv -output_TCR test.csv -output_VJ testVJ.csv -output_log test.log
```
