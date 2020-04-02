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

<img src="https://github.com/jcao89757/TESSA/blob/master/example_data/meta_example_fig.png" width="500">

**Fig.1 |** An example of a TCR meta data matrix in .csv format.

2. A matrix representing the gene expression levels. Columns correspond to cells, rows correspond to genes **(Fig. 2)**. Cells should be in the same order as the cell identifiers in the meta data matrix. Please find the the [.csv example](https://github.com/jcao89757/TESSA/blob/master/example_data/example_exp.csv) for the expression matrix.

<img src="https://github.com/jcao89757/TESSA/blob/master/example_data/exp_example_fig.png" width="500">

**Fig.2 |** An example of an expression data matrix in .csv format.
### Suggested pre-processing workflow
The meta data matrix is suggested to be double-checled and make sure that each element in the column 'contig_id' is unique. The 'cdr3' column allows deplicates, but sequences with any letters, numbers or symbols that do not represent amino acids should be removed. 

For the expression data, the log-transformation is always suggested if the data distribution is right-skewed. Users are encouraged to select their own way to normalize the data. Some useful papers are listed below for you to refer.
1. [Hafemeister, Christoph, and Rahul Satija. "Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression." Genome Biology 20, no. 1 (2019): 1-15.](https://link.springer.com/article/10.1186/s13059-019-1874-1)
2. [Vallejos, Catalina A., Davide Risso, Antonio Scialdone, Sandrine Dudoit, and John C. Marioni. "Normalizing single-cell RNA sequencing data: challenges and opportunities." Nature methods 14, no. 6 (2017): 565.](https://www.nature.com/articles/nmeth.4292)

Moreover, the expression matrix can be replaced with the results from treating it with typical dimension reduction methods, for example, PCA, t-SNE, or UMAP, as long as the format matches the original expression matrix (coordinates or features on the rows, and cell identifiers on the columns).
### Model parameters
The tessa model takes a series of input items listed in the table below.

|Parameters|Description|
|----------|-------|
|tcr|A .csv file contains the meta data matrix described in **Input data**.|
|model|A .h5 file contains a well-trained auto-encoder model used in the BriseisEncoder. Users can find the model in the BriseisEncoder folder.|
|embeding_vectors|A .csv file contains the 'Atchley factors' of all amino acids used in the first step of the BriseisEncoder. Users can find the file in the BriseisEncoder folder.|
|output_TCR|A name of a .csv file that will be created by the model to save the embedded TCRs.|
|output_VJ|(Optional) A name of a .csv file that will be created by the model. If the meta data matrix contains one or two of the optional columns, 'v_gene' and 'j_gene', which denotes the V gene and J gene subgroups (TRBV1-30 and TRBJ1/2) of TRB recombinants, the BriseisEncoder will perform one-hot encoding on those genes and save the encodings in this file for any future uses.|
|output_log|A plain text log file to record any errors or warnings from the BriseisEncoder.|
|exp|A .csv file contains the expression matrix described in **Input data**.|
|output_tessa|A path to save the results data generated by the tessa model.|
|within_sample_networks|A binary value indicating whether the TCR networks are constructed only within TCRs from the same sample/patient (TRUE) or with all the TCRs in the meta data matrix (FALSE). If this is setted to 'TURE', the meta data matrix should include a character column with the name 'sample', which indicates the samples from where the TCRs were generated. The meta data matrix example in **Fig. 1** contains an example 'sample' column.|
|predefined_b|(Optional) A .csv file contains a pre-defined b vector. Please check the paper of tessa for more details about the b vector. If this file is provided, the tessa will not update b in the MCMC iterations. Please find on example for the file [here](https://github.com/jcao89757/TESSA/blob/master/example_data/fixed_b.csv). |
### TCR network construction
The tessa model construct TCR networks with the following command.
```{shell}
python3 Tessa_main.py -tcr ./example_data/example_TCRmeta.csv -model ./BriseisEncoder/TrainedEncoder.h5 -embeding_vectors ./BriseisEncoder/Atchley_factors.csv -output_TCR test.csv -output_VJ testVJ.csv -output_log test.log -exp ./example_data/example_exp.csv -output_tessa /path/to/tessa_results/ -within_sample_networks FALSE
```
After the script finished running, the embedded TCRs can be checked from the saved .csv file. A typical example is shown [here](https://github.com/jcao89757/TESSA/blob/master/example_data/example_TCRembedding.csv) and in **Fig. 3**.

<img src="https://github.com/jcao89757/TESSA/blob/master/example_data/embedding_example_fig.png" width="500">

**Fig.3 |** An example of embedded TCRs in .csv format.

The TCR network result 'tessa_final.RData' is saved in the tessa result folder and can be checked with the following code.
```{r}
load('tessa_final.RData')
m=tessa_results$meta
```
The matrix m has three columns. The 'barcode' column contains the same cell identifiers we used in the expression matrix (the column names) and the meta data matrix (the 'contig_id' column). The 'group_ID's are the TCR sequences of the cells. The column 'cluster_number' contains the TCR sequences of the centered TCRs in the networks. The cells that have the same 'cluster_number' are in the same network. 
## Note 1: Run BriseisEncoder to generate TCR embeddings.
The BriseisEncoder can be used as part of the tessa model, or as a separate tool to generate numerical TCR embeddings for the usage of any other future algorithms, as shown in the code below. All the parameter settings in the following code are the same as we described in **Model parameters**. The embedded TCRs generated by the BriseisEncoder are similar to **Fig. 3**.
```{Shell}
python3 ./BriseisEncoder/BriseisEncoder.py -tcr ./example_data/example_TCRmeta.csv -model ./BriseisEncoder/TrainedEncoder.h5 -embeding_vectors ./BriseisEncoder/Atchley_factors.csv -output_TCR test.csv -output_VJ testVJ.csv -output_log test.log
```
## Note 2: Run tessa with prepared embedded TCRs.
TCR networks can be constructed with the embedded TCRs generated by other algorithms, as shown in the code below. In the **Guided tutorial**, we used our own 30-digits TCR embedding. However, users are free to use any other embeddings of the TCRs with different dimensions, and our software implementation has taken this flexibility in input into consideration. All the parameter settings in the following code are the same as we described in **Model parameters**. Besides, two additional matrices are required in this code. The 'embedding' parameter indicates a .csv file containing the user-generated embedded TCRs, which are in the same format as the embeddings in the [example_TCRembedding.csv](https://github.com/jcao89757/TESSA/blob/master/example_data/example_TCRembedding.csv) (column numbers could be different, and the column names do not matter). The 'meta' parameter indicates a .csv file containing the initial TCR sequences and cell identifiers, which are in the same format as the [meta data matrix](https://github.com/jcao89757/TESSA/blob/master/example_data/example_TCRmeta.csv).
```{shell}
python3 Tessa_main.py -exp ./example_data/example_exp.csv -embedding ./example_data/example_TCRembedding.csv -meta ./example_data/example_TCRmeta.csv -output_tessa /path/to/tessa_results/ -within_sample_networks TRUE -predefined_b ./example_data/fixed_b.csv
```
## Version update
1.0.0: First release. (03-29-2020)



