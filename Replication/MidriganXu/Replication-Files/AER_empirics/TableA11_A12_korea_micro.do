
/**********micro-moments for Korea,  China/Colombia can be defined similarly**************/

clear

clear matrix

set mem 500m

cd "C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\korean_data"

use korea_data_clean_final_wb


**0.  Number of Observations, and Unique Plants

tab year
unique plant

*0.1 Define age/cohort
replace yborn=. if yborn<0|yborn>96
*replace yborn=. if yborn<0|yborn>92

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


*0.2 Start to define alternative group of plants

**based on age
gen cage=1 if age>=0&age<=5
replace cage=2 if age>5&age<=10
replace cage=3 if age>10&age~=.

** 0.3 construct plant level tfp with index number
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

gen mal2=mal*(0.66/0.52)
gen tfp=(lnvp-mal2*0.85*lnl-(1-mal2)*0.85*lnk)/(1-0.85)
sum tfp, d
replace tfp=. if tfp>r(p99)|tfp<r(p1)

xtset plant year

**generate De
gen De=tfp-L.tfp

*take the mean of tfp as z
egen zp=mean(tfp), by (plant)

egen nobs=count(year), by (plant)
gen e=tfp-zp

***define which group we want to use
gen group=cage
tab group, gen (dgroup)

gen vpk=lnvp-lnk
sum vpk, d
replace vpk=. if vpk>r(p99)|vpk<r(p1)

********************************var/growth rate for alternative variables**************
foreach var of varlist vpk lnv lnk tfp{

**level
sum `var', d
matrix s`var'=[r(sd)]
matrix m`var'=[r(mean)]
tabstat `var', by (group) stat(sd) save
matrix s`var'a=[r(Stat1), r(Stat2), r(Stat3)]
tabstat `var', by (group) stat(mean) save
matrix m`var'a=[r(Stat1), r(Stat2), r(Stat3)]

**differencing
gen d`var'=`var'-L.`var'
sum d`var',d
matrix sd`var'=r(sd)
matrix md`var'=r(mean)
tabstat d`var', by (group) stat (sd) save
matrix sd`var'a=[r(Stat1), r(Stat2),r(Stat3)]
tabstat d`var', by (group) stat (mean) save
matrix md`var'a=[r(Stat1), r(Stat2), r(Stat3)]
}



*************************************elast. D(Y or K) on De  ******************
gen D1e=De*dgroup2
gen D2e=De*dgroup3

gen d1lnv=dlnv*dgroup2
gen d2lnv=dlnv*dgroup3

foreach var of varlist dlnv dlnk dvpk {

areg  `var' De, absorb(indy)
matrix re`var'=[_coef[De]]

areg `var' De D1e D2e, absorb(indy)
matrix re`var'a=[_coef[De],_coef[De]+_coef[D1e],_coef[De]+_coef[D2e]]

areg `var' dlnv, absorb(indy)
matrix rv`var'=[_coef[dlnv]]

areg `var' dlnv d1lnv d2lnv, absorb(indy)
matrix rv`var'a=[_coef[dlnv],_coef[dlnv]+_coef[d1lnv],_coef[dlnv]+_coef[d2lnv]]

}

matrix datam=[sdlnv'\sdlnva'\sdlnk'\sdlnka'\redlnv'\redlnk'\rvdlnk'\mdlnva'\mdlnka'\mvpka'\svpk'\svpka']

************************distribution of growth rate by age***************
gen dlnvp=lnvp-L.lnvp
tabstat dlnvp, by (group) stat (p25 p50 p75 sd)
