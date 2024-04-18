clear

clear matrix

set mem 500m

cd "C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\colombia_data"

use col_data_clean_final_wb


**0.  Number of Observations, and Unique Plants
replace sic=floor(sic/10) if sic>10000
gen inid2=floor(sic/10)

egen nyear=count(year), by (plant)
drop if nyear>11
sort plant year
by plant: gen error_y=1 if year==year[_n-1]
drop if error_y==1

/*egen minyear=min(year), by (plant)
egen maxyear=max(year), by (plant)*/

*0.1 Define age/cohort
replace yborn=. if yborn<0|yborn>91

*make sure birth year same for the same firm
egen by_nobs=count(year), by (plant yborn)
egen by_mobs=max(by_nobs), by (plant)
gen  dmost=1 if by_nobs==by_mobs
*define cohort
egen cohort=min(yborn*dmost), by (plant)
drop by_nobs by_mobs dmost 

gen age=year-cohort
**check miss reporting
replace age=. if age<0
**
replace age=. if age==0&cohort>minyear
count if age==.

keep if year>=85&year<=90
tab year
unique plant

/**construct plant level tfp with index number**/
**
gen lnl2=lnl^2
gen lnl3=lnl^3
gen lnk2=lnk^2
gen lnk3=lnk^3
gen lnlk=lnl*lnk
gen lnl2k=lnl2*lnk
gen lnlk2=lnl*lnk2
** purified lnv
egen indy=group(inid2 year)
areg lnv lnl lnk lnl2 lnk2 lnl3 lnk3 lnlk lnl2k lnlk2, absorb(indy)
predict lnvp, xb
**input elas. of labor
gen lnal=lnl-lnv
sum lnal, d
replace lnal=. if lnal>r(p95)|lnal<r(p5)
egen mlnal=mean(lnal), by (inid2)
gen mal=exp(mlnal)/0.85

gen tfp1=lnv-mal*lnl-0.85*(1-mal)*lnk
gen tfpp1=lnvp-mal*lnl-0.85*(1-mal)*lnk

gen mal2=mal*(0.66/0.51)
gen tfp=(lnvp-mal2*0.85*lnl-(1-mal2)*0.85*lnk)/(1-0.85)
sum tfp, d
replace tfp=. if tfp>r(p99)|tfp<r(p1)

tabstat mal2, by (inid2)
