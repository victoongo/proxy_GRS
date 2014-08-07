set more off
set matsize 11000


***** Step A, preping the marker info files to get the rs snp names
*ARIC
cd D:\data_work\aric\geno
insheet using "D:\data_work\aric\geno\marker_info\GenomeWideSNP_6.na27.annot.csv", names clear
rename dbsnprsid rs
sort rs
save GenomeWideSNP_6.na27.annot.dta, replace

	***** Step A1, extract probesetid and rs for later conversion
	use GenomeWideSNP_6.na27.annot.dta, clear
	* 2. KEEP ONLY RS NAMES THAT WE WANT - DROP DUPLICATES (REQUIRES EXTRA STEP TO PICK THE BEST)
	keep probesetid rs
	bysort rs: keep if _n==1 //(2987 observations deleted)
	save crossref_list, replace
	* 3. EXTRACT REMAINING PROBESET ID TO TEXT FILE
	outfile probesetid using "probesetid_list.raw", noquote replace
	* 4. EXTRACT SAME SNPS WITH RS NAMES TO ANOTHER TEXT FILE
	outfile probesetid rs using "probesetid_rs_list.raw", noquote replace

*CHS
cd D:\data_work\chs\geno
/*
insheet using "D:\data_work\chs\geno\marker_info\HumanCNV370v1-Duo.csv", comma clear
gen nhash=strmatch(v1, "#*")
drop if nhash==1
drop nhash
do "D:\victor\Dropbox\P_dbGaP\first_row_as_varname.do"
sort rs
save HumanCNV370v1-Duo.dta, replace
*/
insheet using "D:\data_work\chs\geno\marker_info\HumanOmni1-Quad_v1-0_B.csv", clear
drop if _n<=7
do "D:\victor\Dropbox\P_dbGaP\first_row_as_varname.do"
rename name rs
sort rs
save HumanOmni1-Quad_v1-0_B.dta, replace

insheet using "D:\data_work\chs\geno\marker_info\HumanOmni1-Quad_v1-0_H.csv", clear
drop if _n<=7
do "D:\victor\Dropbox\P_dbGaP\first_row_as_varname.do"
rename name rs
sort rs
save HumanOmni1-Quad_v1-0_H.dta, replace

*MESA
cd D:\data_work\mesa\geno
insheet using "D:\data_work\mesa\geno\marker_info\GenomeWideSNP_6.na24.annot.csv", clear
drop if _n<=15
do "D:\victor\Dropbox\P_dbGaP\first_row_as_varname.do"
rename dbsnp_rs_id rs
sort rs
save GenomeWideSNP_6.na24.annot.dta, replace

	***** Step A1, extract probesetid and rs for later conversion
	use GenomeWideSNP_6.na24.annot.dta, clear
	keep probe_set_id rs
	bysort rs: keep if _n==1 //(10278 observations deleted)
	save crossref_list, replace
	outfile probe_set_id using "probesetid_list.raw", noquote replace
	outfile probe_set_id rs using "probesetid_rs_list.raw", noquote replace
	

*SAGE
cd D:\data_work\sage\geno
insheet using "D:\data_work\sage\geno\marker_info\SNP_annotation.csv", names clear
rename rsid rs
sort rs
save SNP_annotation.dta, replace


***** Step B, merge npu and gru 
* CHS
cd "D:\data_work\chs\geno"
!plink --bfile original\c1_CHS_v3 --bmerge original\c2_CHS_v3.bed original\c2_CHS_v3.bim original\c2_CHS_v3.fam --make-bed --out CHS_v3

* MESA
cd "D:\data_work\mesa\geno"
!plink --bfile original\SHARE_MESA_c1 --bmerge original\SHARE_MESA_c2.bed original\SHARE_MESA_c2.bim original\SHARE_MESA_c2.fam --make-bed --out share_mesa


***** Step C, convert probesetid to RS
* ARIC
cd D:\data_work\aric\geno
* 5. EXTRACT BINARY DATA BASED ON PROBESET ID
!plink2 --bfile ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out --extract probesetid_list.raw --make-bed --out probesetid_list
* 6. USE UPDATE NAMES TO CHANGE THE NEW DATA'S PROSETID NAMES TO RS NAMES
!plink2 --bfile probesetid_list --update-map probesetid_rs_list.raw --update-name --make-bed --out ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_rs

* MESA
cd D:\data_work\mesa\geno
* 5. EXTRACT BINARY DATA BASED ON PROBESET ID
!plink2 --bfile share_mesa --extract probesetid_list.raw --make-bed --out probesetid_list
* 6. USE UPDATE NAMES TO CHANGE THE NEW DATA'S PROSETID NAMES TO RS NAMES
!plink2 --bfile probesetid_list --update-map probesetid_rs_list.raw --update-name --make-bed --out share_mesa_rs

***** Step D, convert ped to binary 
* Dunedin
cd "D:\data_work\dunedin\geno"
!plink2 --file Dunedin_SNP_arrays_IDcorrected --make-bed --out Dunedin_SNP_arrays_IDcorrected



***** extract
*ARIC
cd D:\data_work\aric\geno
use GenomeWideSNP_6.na27.annot.dta, clear
merge m:1 rs using "D:\victor\data_live\extract_lst.dta", nogen keep(3)
keep probesetid rs
bysort rs: gen n=_n
save extract_lst_aric, replace
outfile probesetid using "extract_lst_aric.txt", noquote replace
*!plink2 --bfile extracted --geno 0.05 --maf 0.05 --hwe 0.001 --recode --out extracted_qc 
!plink2 --bfile ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out --extract extract_lst_aric.txt --recode --out extracted
!plink2 --file extracted --recode-lgen --out extracted_l 
*!plink2 --file extracted --geno 0.05 --maf 0.05 --hwe 0.001 --recode --out extracted_qc 

insheet using extracted_l.lgen, delim(" ") clear
rename (v1 v2 v3 v5 v6) (fid iid probesetid a1 a2)
drop v4
sort probesetid
save extracted_l, replace

use extracted_l, clear
merge m:1 probesetid using extract_lst_aric, nogen
sort rs
merge m:m rs using "D:\victor\data_live\extract_lst_info.dta", nogen keep(1 3)
gen proxy=1 if rs~=snp
save extracted_final, replace



*MESA
cd D:\data_work\mesa\geno
use GenomeWideSNP_6.na24.annot.dta, clear
merge m:1 rs using "D:\victor\data_live\extract_lst.dta", nogen keep(3)
keep probe_set_id rs
bysort rs: gen n=_n
save extract_lst_mesa, replace
outfile probe_set_id using "extract_lst_mesa.txt", noquote replace
!plink2 --bfile share_mesa --extract extract_lst_mesa.txt --recode --out extracted
!plink2 --file extracted --recode-lgen --out extracted_l 
*!plink2 --file extracted --geno 0.05 --maf 0.05 --hwe 0.001 --recode-lgen --out extracted_qc 

insheet using extracted_l.lgen, delim(" ") clear
rename (v1 v2 v3 v5 v6) (fid iid probe_set_id a1 a2)
drop v4
sort probe_set_id
save extracted_l, replace

use extracted_l, clear
merge m:1 probe_set_id using extract_lst_mesa, nogen
sort rs
merge m:m rs using "D:\victor\data_live\extract_lst_info.dta", nogen keep(1 3)
gen proxy=1 if rs~=snp
save extracted_final, replace

*SAGE
cd D:\data_work\sage\geno
use SNP_annotation.dta, clear
merge m:1 rs using "D:\victor\data_live\extract_lst.dta", nogen keep(3)
keep rs
bysort rs: gen n=_n
save extract_lst_sage, replace
outfile rs using "extract_lst_sage.txt", noquote replace
!plink2 --bfile GENEVA_SAGE_HR_GENO_FINAL_ZEROED --extract extract_lst_sage.txt --recode --out extracted
!plink2 --file extracted --recode-lgen --out extracted_l 
*!plink2 --file extracted --geno 0.05 --maf 0.05 --hwe 0.001 --recode --out extracted_qc 

insheet using extracted_l.lgen, delim(" ") clear
rename (v1 v2 v3 v5 v6) (fid iid rs a1 a2)
drop v4
sort rs
save extracted_l, replace

use extracted_l, clear
merge m:1 rs using extract_lst_sage, nogen
sort rs
merge m:m rs using "D:\victor\data_live\extract_lst_info.dta", nogen keep(1 3)
gen proxy=1 if rs~=snp
save extracted_final, replace

***** Quality Control


*****

