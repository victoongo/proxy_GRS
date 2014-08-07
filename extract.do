



use mesalst$snplst, clear
merge 1:m rs using markerinfo, nogen keep(1 3)
save mesalst$snplsta, replace
use mesalst$snplsta, clear
bysort loc_snp_id: keep if _n==1
save mesalst1b, replace
keep loc_snp_id
bysort loc_snp_id: keep if _n==1
outsheet using mesalst1.txt, nonames noquote replace

!plink --bfile npu\SHARE_MESA_c2 --missing --out mesac1
!plink --bfile cohort\SHARE_MESA_c1 --missing --out mesac2

!plink --bfile npu\SHARE_MESA_c2 --extract mesalst1.txt --recode --out mesac2lst1
!plink --bfile cohort\SHARE_MESA_c1 --extract mesalst1.txt --recode --out mesac1lst1

!plink --file mesac1lst1 --merge mesac2lst1.ped mesac2lst1.map --recode --out mesaclst1

!plink --file mesaclst1 --recodeA --out mesaclst1

!plink --file mesaclst1 --missing --out mesaclst1
!plink --file mesaclst1 --freq --out mesaclst1

insheet using "mesaclst1.map", tab clear
keep v2
gen n=_n
rename v2 loc_snp_id
merge 1:1 loc_snp_id using mesalst1b, nogen
sort n
keep rs
local snplst 
forvalues n= 1/`c(N)' {
	local snplst `snplst' `= rs[`n']'_a1 `= rs[`n']'_a2 
}
di "`snplst'"
global snplst "`snplst'"
di "$snplst"
*insheet fid iid mo fa sex pheno snplst $snplst using "mesaclst1.ped", delim(" ") clear
* insheet doesn't allow for more than 256 charactors. 

insheet using "mesaclst1.ped", delim(" ") clear
unab vlst: v*
global vlst `vlst'
di "$vlst"
rename ($vlst) (fid iid mo fa sex pheno $snplst)
save mesaclst1_f, replace
