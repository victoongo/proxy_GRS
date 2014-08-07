set more off

***** geno data
* 
do "proxysearch_readexcel.do"
***** Manual step, upload the input_list.raw file to the snap interface, and select according to the setttings in SNAPsettings.PNG 
***** Manual step, move download from SNAP to victor/datalive. Name should be "SNAPresults.txt" 
do "proxysearch.do"
do "geno.do"
do "extract.do"

***** pheno data
do "pheno.do"
