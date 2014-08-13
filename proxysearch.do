cd "D:\victor\data_live"
set more off

***** run proxysearch_readexcel first each time you run this to get the pheno list into a global *****


***** Step B: clean the SNAP results and prep for SNP extraction from downloaded geno data
local ref /*hm22 hm3*/ 1000g 
*local refcount : word count `ref' //Use if multiple SNAP files (from 1000g, HM2, HM3)
local refcount : list sizeof local(ref) //Use if only 1 SNAP file (from 1000g)
di `refcount'
if `refcount'==1 {
	insheet using "SNAPResults.txt", clear
	gen ref="`ref'"
	drop if proxy=="WARNING"
	save SNAPResults_all, replace
}
else {
	foreach x of local ref {
		insheet using "SNAPResults_`x'.txt", clear
		gen ref="`x'"
		save SNAPResults_`x', replace
	}
	forval i=1/`refcount' {
		if `i'==1 {
			use SNAPResults_`: word `i' of `ref'', clear
		}
		else {
			append using SNAPResults_`: word `i' of `ref''
		}
	}
}
save SNAPResults_all, replace

use SNAPResults_all, clear
bysort snp proxy ref: drop if _n>1
destring *, replace
quietly: ta ref
local n_ref=r(r)
if `n_ref'>1 {
	/*
	gen b_=1
	local chgvlst distance rsquared 
	*local chgvlst_new
	foreach x of local chgvlst  {
		rename `x' `x'_
		*local chgvlst_new `chgvlst_new' `x'_
	}
	reshape wide b_ distance rsquared , i(snp proxy) j(ref, string)
	*/
}
/*
gen b_=1
reshape wide b_, i(snp proxy) j(ref, string)
* reshape wide produced this error, which means different sources produced slightly different results. 
* this means 
distance not constant within snp proxy
rsquared not constant within snp proxy
dprime not constant within snp proxy
*/
/*
local arrays OQ A6 IM
foreach x of local arrays {
	gen `x'=1 if strpos(arrays, "`x'")
	*keep if strmatch(arrays,"*`x'*")
}
* see arrays names here: http://www.broadinstitute.org/mpg/snap/ldsearch.php
*/
rename snp rs
sort rs
gen n=1
merge m:1 rs n using snplists_w
foreach p of global phenolst {
	di "`p'"
	gen prx_ra_`p'=""
	*** with flag
	replace prx_ra_`p'=minor if (flag_maf45_`p'==1 | flag_atcg_`p'==1) & ra_`p'==minor_`p'
	replace prx_ra_`p'=major if (flag_maf45_`p'==1 | flag_atcg_`p'==1) & ra_`p'==major_`p'
	*** without flag
	gen strand_r=.
	replace strand_r=1 if (ra_r_`p'==minor) & rs==proxy
	replace strand_r=2 if (ra_r_`p'==major) & rs==proxy
	bysort rs: egen strand_rmean=mean(strand_r)
	ta strand_rmean
	replace prx_ra_`p'=minor if (flag_maf45_`p'~=1 & flag_atcg_`p'~=1) & strand_rmean==1
	replace prx_ra_`p'=major if (flag_maf45_`p'~=1 & flag_atcg_`p'~=1) & strand_rmean==2
	replace prx_ra_`p'=ra_`p' if (flag_maf45_`p'~=1 & flag_atcg_`p'~=1) & strand_rmean==.
	drop strand_r strand_rmean
}
sort rs distance
drop arrays
*** // need to pick the best for each phenotype? 

* SAVE TEXT FILE WITH RS NAMES ONLY
outfile proxy using "proxy_list.raw", noquote replace
* SAVE STATA DATA WITH RS NAMES AND RISK ALLELE INFORMATION (DATA_PROXY)
keep rs proxy distance rsquared es_* prx_ra_*
drop if proxy==""
rename (rs proxy) (origin_rs rs)
save prx_ra, replace
* SAVE A TRANSPOSED VERSION // no need for this if use method 2 below


***** Step C: prep the downloaded geno data for extraction
cd "D:\victor\data_live"
local studies dunedin aric mesa  //sage
local aric ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_rs
local mesa SHARE_MESA_rs
*local sage GENEVA_SAGE_HR_GENO_FINAL_ZEROED_qc
local dunedin Dunedin_SNP_arrays_IDcorrected
*local chs CHS_v3

*global phenolst HDL // for testing only

pause on
foreach x of local studies {
	di "`x'"
	foreach y of local `x' {
		di "`y'"
		*** method 1
		/*
		insheet using "D:\data_work/`x'\geno/`y'.bim", tab clear // NOW ALL IN RS NAMES
		*MERGE WITH THE PROXY LIST //KEEP THE SNP NAMES THAT ARE IN BOTH LIST
		merge 1:m rs using "D:\victor\data_live\prx_rs.dta", keep(3) 
		*SAVE AS TEXT FILE
		keep rs
		outfile proxy using "proxy_list.raw", noquote replace

		*USE PLINK TO EXTRACT DATA BASED ON THE RS NAMES (ALSO CONVERTED TO PED/MAP FILES)
		!plink2 --bfile share_mesa --extract probesetid_list.raw --recode --out probesetid_list
		*STATA READ MAP FILES TO CREATE THE DO FILE SO NAME THE DATA FROM PED
		*USE STATA TO READ THE PED FILES AND USE THE ABOVE DO FILE TO NAME THE COLUMNS
		
		*MERGE WITH DATA_PROXY(TRANSPOSED) AND CALCULATE RISK SCORE FOR EACH SNP
		!plink2 --bfile share_mesa --extract probesetid_list.raw --recode --out probesetid_list
		*/
		*** method 2
		local c=0
		foreach p of global phenolst {
			di "`p'"
			* open the prx_ra file and keep only unique ones
			use prx_ra, clear
			
			keep origin_rs rs distance rsquared es_`p' prx_ra_`p'
			
			order rs prx_ra_`p' es_`p' 
			drop if rs=="" | es_`p'==. | prx_ra_`p'==""
			
			merge m:1 rs using D:\data_work/`x'\geno/rs_list, nogen keep (3)
			if _N>0 {
				
				* create the proxy snp list for extraction // if you want to keep all the proxy snps, use the the following line to extract all proxies
				*outfile rs prx_ra_`p' using "D:\data_work/`x'\geno/prx_ra_`p'.raw", noquote replace 
				
				** begin: select one proxy for each origial snp
				gen prx_f=0
				replace prx_f=1 if distance==0
				bysort origin_rs: egen has_origin=max(prx_f)
				drop if prx_f==0 & has_origin==1
				bysort origin_rs: egen high_rsq=max(rsquared) if distance>0
				keep if prx_f==1 | rsquared==high_rsq
				bysort origin_rs: egen shortdist_rsq=min(distance) if rsquared==high_rsq & distance>0
				keep if prx_f==1 | distance==shortdist_rsq 
				** end
				
				if _N>0 {
					save D:\data_work/`x'\geno/proxy_`p', replace 
					*bysort rs: keep if _n==1
					* output rs and ra and es as text file
					keep rs prx_ra_`p' es_`p'
					order rs prx_ra_`p' es_`p'
					outfile using "D:\data_work/`x'\geno/prx_es_ra_`p'.raw", noquote replace
					outfile rs using "D:\data_work/`x'\geno/prx_`p'.raw", noquote replace
					
					** begin: create the risk allele count file for stata
					* create the proxy snp list for extraction // the following line only keep one proxy for each original rs. 
					* // if you want to keep all the proxy snps, comment out the following line and use the line before the "select one proxy for each snp" code to extract all proxies
					outfile rs prx_ra_`p' using "D:\data_work/`x'\geno/prx_ra_`p'.raw", noquote replace
					* use plink to extract the snps
					!plink2 --bfile D:\data_work/`x'\geno/`y' --extract "D:\data_work/`x'\geno/prx_`p'.raw" --recode --out D:\data_work/`x'\geno/prx_`p'
					* plink to calculate ref allele count
					!plink --file D:\data_work/`x'\geno/prx_`p' --recodeA --recode-allele "D:\data_work/`x'\geno/prx_ra_`p'.raw" --out D:\data_work/`x'\geno/prx_ac_`p' // plink2 doesn't work here and --reference-allele is different
					* read the risk allele count into stata and recode the missing to .
					insheet using "D:\data_work/`x'\geno/prx_ac_`p'.raw", delimiter(" ") clear
					foreach snpname of varlist rs* {
						destring `snpname', replace ignore(NA)
					}
					save "D:\data_work/`x'\geno/prx_ac_`p'.dta", replace
					** end
					
					* plink to calculate GRS with --score (rs ra es)
					*!plink2 --bfile D:\data_work/`x'\geno/prx_`p' --score D:\data_work/`x'\geno/prx_es_ra_`p'.raw --out D:\data_work/`x'\geno/grs_`p' // this will give a different genotyping rate in the log file
					!plink2 --bfile D:\data_work/`x'\geno/`y' --score D:\data_work/`x'\geno/prx_es_ra_`p'.raw --out D:\data_work/`x'\geno/grs_`p'
					
					local profile : dir "D:\data_work/`x'\geno/" files "grs_`p'.profile", respectcase
					if `"`profile'"'~="" {
						insheet using D:\data_work/`x'\geno/grs_`p'.profile, clear
						drop if _n==1
						gen pheno_name="`p'"
						if `c++'>0 append using D:\data_work/`x'\geno/grs_all, force
						save D:\data_work/`x'\geno/grs_all, replace
					}
				}
			}
		}
		*des
		*sum
		*insheet using "D:\data_work/`x'\geno/`y'.fam", delimiter(" ") clear
	}
	local grs : dir "D:\data_work/`x'\geno/" files "grs_all.dta", respectcase
	if `"`grs'"'~="" {
		use D:\data_work/`x'\geno/grs_all, clear
		forvalues c=1/10 {
			replace v1=subinstr(v1,"  "," ",.)
		}
		replace v1=trim(v1)
		split v1, p(" ")
		drop v1 v13
		rename (v11 v12 v14 v15 v16) (fid iid cnt cnt2 score)
		order fid iid cnt cnt2 score pheno_name
		save D:\data_work/`x'\geno/grs_all, replace
	}
}
* verify data and choose zeroed out version???
* may need to combine data such as GPU and NPU
* data quality check


