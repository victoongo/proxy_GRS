unab varnames: *
foreach y of local varnames {
		quietly: ta `y'
		local n=r(r)
		*di `n'
		if `n'>0 {
			replace `y'=lower(`y') if _n==1
			quietly replace `y'=subinstr(`y', " ", "_", .) if _n==1 
			quietly replace `y'=subinstr(`y', ".", "_", .) if _n==1 
			quietly replace `y'=subinstr(`y', "-", "_", .) if _n==1 
			quietly replace `y'=subinstr(`y', "/", "_", .) if _n==1 
			quietly replace `y'=subinstr(`y', "#", "", .) if _n==1 
			if length(`y'[1])>32 {
				replace `y'=abbrev(`y', 32) if _n==1 
			}		
			rename `y' `= `y'[1]'
		}
		else drop `y' // drop variable if the variable has no valid values. 
}
drop if _n==1
destring *, replace
