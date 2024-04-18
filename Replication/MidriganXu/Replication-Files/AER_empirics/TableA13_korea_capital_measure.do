

*********CODE TO PRODUCE KOREAN  DATA********

****note: benchmark use machinery_equipment_building

***drop observations only necessary - when value-added, labor, or capital measures<=0, if drop, drop the whole plant

clear

clear matrix

set mem 1g

cd "C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\korean_data"

use "C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\korean_data\raw data processing\korea_data.dta"

sort busstr year
egen plant=group(busstr)

sort year inid2

keep if year<=96

/*construct the capital stock and investment flow using perpetual inventory method*/

***without building
/*gen ecap_b=(mchby+carby)/iprice   /**book value BOY for comparison**/
gen ecap_p=(mchby+carby)/iprice

gen einv=(mchpch+carpch-mcsell-casell)/iprice
gen eret=(mchdep+cardep)/iprice*/

***with building
gen ecap_b=(bldgby+mchby+carby)/iprice   /**book value BOY for comparison**/
gen ecap_p=(bldgby+mchby+carby)/iprice

gen einv=(bdpch+mchpch+carpch-bdsell-mcsell-casell)/iprice
gen eret=(bddep+mchdep+cardep)/iprice

sort busstr year
by busstr: egen minyear=min(year)
by busstr: egen maxyear=max(year)

by busstr: replace ecap_p=ecap_p[_n-1]+einv[_n-1]-eret[_n-1] if year~=minyear

tabstat ecap_p, by (year) stat (p1 p5 p25 p50 p75 p95)
tabstat ecap_b, by (year) stat (p1 p5 p25 p50 p75 p95)

corr ecap_p ecap_b

/*construct material, labor, capital, and output measures, all deflated by CPI*/

*total real revenue
gen vq=produt
gen q=produt/cpi

*material input  and energy input
gen vm= (mtrl+maint)+subout+(fuel)+eltry+wtr
gen ve=(fuel)+eltry+wtr
gen qm=vm/cpi
gen qe=ve/cpi

*value-added
gen am=vm/produt
replace am=. if am>1|am<0
egen mam=mean(am), by (inid2 year)

*labor input 
gen ql=(twage+wfare)/cpi

*capital and investment
gen qk=ecap_p+rent/cpi/uk
gen qi=einv

*gen qk=(fuel+eltry)/cpi    /***proxy used by energy***/

/*construct firm level log(inputs) and log(output)*/
gen lnq= ln(q)
gen lnl = ln(ql)
gen lnk = ln(qk) 
gen lnm= ln(qm)

*log value-added
gen vv=produt-vm
gen v=vv/cpi
gen lnv=ln(v)

/**drop outliers**/
*drop ones with negative capital stocks 
gen error_cap=1 if qk<=0
replace error_cap=0 if error_cap==.

*drop ones with negative value-added, total labor costs
gen error_rev=1 if v<=0|ql<=0
replace error_rev=0 if error_rev==.

gen dumdrop=error_cap+error_rev     
egen flagdrop=sum(dumdrop), by (busstr)
replace flagdrop=1 if flagdrop>1

**now calculate share of output by plants that dropped each year
table year flagdrop, c(sum produt)

drop if flagdrop==1
unique plant

save korea_data_clean_final_wb.dta, replace
