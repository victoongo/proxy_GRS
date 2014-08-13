set more off

***** prep geno data

cd "C:\Users\dan\Dropbox\Research\P_dbGaP\do"
do "geno.do"

***** prep SNP request data
 
cd "C:\Users\dan\Dropbox\Research\P_dbGaP\do"
do "proxysearch_readexcel.do"

***** Manual step, upload the input_list.raw file to the snap interface, and select according to the setttings in SNAPsettings.PNG 
***** Manual step, move download from SNAP to victor/datalive. Name should be "SNAPresults.txt" 

cd "C:\Users\dan\Dropbox\Research\P_dbGaP\do"
do "proxysearch.do"




***** pheno data
do "pheno.do"
