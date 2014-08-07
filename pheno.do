set more off
set matsize 11000

*cd "D:\Projects\test"

***** Step A, convert CSV to DTA, with the correct var name and label
local proj mesa/pheno/gru mesa/pheno/npu sage/pheno aric/pheno chs/pheno/gru chs/pheno/npu
foreach p of local proj {
	cd "D:\data_work/`p'"
	local rlst : dir . files "*.txt"
	di `"`rlst'"'
	* option 2: directly rename and label each var. the problem of this approch is the destirng part is very inefficient. 
	foreach x in `rlst' {
		local dlst : dir . files "`x'.dta"
		di `"`dlst'"'
		if `"`dlst'"'=="" { // comment out check for existence
			di "does not exist!!!"
			insheet using `x', clear
			gen nhash=strmatch(v1, "#*")
			replace nhash=2 if strmatch(v1, "##*")==1
			drop if nhash==1
			drop nhash
			unab varnames: *
			foreach y of local varnames {
				quietly tostring `y', replace
				
				*macro dir
				*macro list
				if (`y'[1]=="." | `y'[1]=="") & (`y'[2]=="." | `y'[2]=="") {
					*drop `y' 
				}
				else {
					lab var `y' "`= `y'[1]'"
					quietly replace `y'=subinstr(`y', " ", "_", .) if _n==2 // this changes the space to underscore in all var names.
					quietly replace `y'=subinstr(`y', ".", "_", .) if _n==2 // this changes the . to underscore in all var names.
					rename `y' `= `y'[2]'
				}
			}
			drop if _n<3
			*quietly destring *, replace
			capture destring dbGaP_SubjID, replace
			capture destring dbGaP_SampID, replace
			save `x'.dta, replace
		} // comment out check for existence
	}
}




* option 1: create two data sets with one for data and the other for var name and label to create the program to rename and label 
/*
	insheet using phs000092.v1.pht000637.v1.p1.AlcoholDepAdd_Pedigree.MULTI.txt, clear
	gen nhash=strmatch(v1, "#*")
	replace nhash=2 if strmatch(v1, "##*")==1
	drop if nhash==1
	tempfile original
	save `original', replace
	drop if nhash==2
	drop if _n==1
	tempfile originaldta
	save `originaldta', replace 
	use `orignal', clear
	drop if _n>2
	transpose
	rename
*/


***** Step B, merge GRU and NPU 
local proj /*mesa */chs
foreach p of local proj {
	cd "D:\data_work/`p'/pheno\npu"
	local rlst : dir . files "*.dta"
	di `"`rlst'"'
	*di `rlst'
	foreach x of local rlst {
		di "`x'"
		if strpos("`x'",".c2.")>0 & strpos("`x'","-npu")>0 {
			local b=subinstr("`x'",".c2.",".c1.",.)
			di "`b'" " **b1"
			local a=subinstr("`b'","-npu","-all",.)
			local b=subinstr("`b'","-npu","",.)
			di "`b'" " **b2"
			*di `b'
			local dlst : dir "D:\data_work/`p'\pheno" files "`b'"
			di `"`dlst'"' "!!!"
			if `"`dlst'"'=="" {
				di "does not exist!!!"
			}
			else {
				local alst : dir "D:\data_work/`p'\pheno" files "`a'"
				di `"`alst'"' "!!!"
				if `"`alst'"'=="" { // comment out check for existence
					di "yeah!!!"
					use `x', clear
					*local id dbGaP_SubjID dbGaP_SampID
					*foreach idi of local id {}
					capture confirm variable dbGaP_SubjID
					if !_rc {
						sort dbGaP_SubjID
						merge m:m dbGaP_SubjID using "D:\data_work/`p'\pheno\gru/`b'" //, force 
						quietly destring *, replace
						save "D:\data_work/`p'\pheno/`b'", replace
					}
					else {
						capture confirm variable dbGaP_SampID
						if !_rc {
							sort dbGaP_SampID
							merge m:m dbGaP_SampID using "D:\data_work/`p'\pheno\gru/`b'" //, force 
							quietly destring *, replace
							save "D:\data_work/`p'\pheno/`b'", replace
						}
					}
				} // comment out check for existence
			} 
		}
	}
}
