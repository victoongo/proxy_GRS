cd "D:\victor\data_live"
set more off


***** Step A: prep snp lists for SNAP by getting the snps needed from an excel file 
import excel "D:\victor\Dropbox\P_dbGaP\SNP Lists.xlsx", sheet("list") clear
do "D:\victor\Dropbox\P_dbGaP\do\first_row_as_varname.do"
unab allvar: *
local idvar pheno rs
local varlst: list allvar-idvar
foreach x of local allvar {
	local type: type `x'
	if substr("`type'",1,3)=="str" {
		replace `x'=trim(`x')
	}
}
drop if rs==""
/*
replace allele=subinstr(allele, " ", "", .)
replace allele=subinstr(allele, "/", "", .)
tostring oa ra_r, replace

gen a1=substr(allele,1,1)
gen a2=substr(allele,3,3)
replace oa=a1 if (oa=="" | oa==".") & ra==a2
replace oa=a2 if (oa=="" | oa==".") & ra==a1
*/
gen ra_r=""
replace ra_r="A" if ra=="T"
replace ra_r="T" if ra=="A"
replace ra_r="C" if ra=="G"
replace ra_r="G" if ra=="C"
gen flag_atcg=1 if (ra=="C" & oa=="G") | (ra=="G" & oa=="C") | (ra=="A" & oa=="T") | (ra=="T" & oa=="A")
replace pheno=subinstr(pheno," ","_",.)
replace pheno=subinstr(pheno,"'","_",.)
replace pheno=subinstr(pheno,"(","_",.)
replace pheno=subinstr(pheno,")","_",.)
replace pheno=subinstr(pheno,"/","__",.)
replace pheno=upper(pheno)
replace pheno=substr(pheno,1,18)
replace raf="" if raf=="NA"
destring raf, replace
gen major=ra
replace major=oa if raf<=0.5
gen minor=ra
replace minor=oa if major==ra
gen maf=cond(raf<=0.5,raf,cond(raf>0.5 & raf~=.,1-raf,.))
gen flag_maf45=1 if maf>0.45 & maf~=.
*keep rs pheno
*sort rs
*save snplists_l, replace
unab allvar: *
local varlst: list allvar-idvar
local vlst
foreach x of local varlst {
	rename `x' `x'_
	local vlst `vlst' `x'_
}
gen b_=1
quietly: levels pheno, local(phenolst)
global phenolst `"`phenolst'"'
di `"$phenolst"'

foreach p of global phenolst {
	di "`p'"
}
bysort rs: gen n=_n
reshape wide b_ `vlst', i(rs n) j(pheno, string)
rename b_* *
sort rs
save snplists_w, replace
outfile rs using "input_list.raw", noquote replace

