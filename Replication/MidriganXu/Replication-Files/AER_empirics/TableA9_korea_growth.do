
/***China and Colombia patterns are constructed in the same way***/

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

*take the mean of tfp as z
egen zp=mean(tfp), by (plant)

egen nobs=count(year), by (plant)
gen e=tfp-zp

gen vpk=lnvp-lnk
sum vpk, d
replace vpk=. if vpk>r(p99)|vpk<r(p1)

xtset plant year

areg vpk, absorb(indy)
predict shifter, d
gen vpkk=vpk-shifter

***
xtset plant year
gen dlnvp=lnvp-L.lnvp
sum dlnvp, d
gen dv75=r(p90)
gen dv25=r(p10)

tabstat vpk if dlnvp>dv75&dlnv~=., stat(var mean)
tabstat vpk if dlnvp<=dv25&dlnv~=., stat(var mean)

gen dlnk=lnk-L.lnk
sum dlnk, d
gen dk75=r(p90)
gen dk25=r(p10)
tabstat vpk if dlnk>=dk75&dlnk~=., stat(var mean)
tabstat vpk if dlnk<=dk25&dlnk~=., stat(var mean)

gen dlne=tfp-L.tfp
sum dlne, d
gen de75=r(p90)
gen de25=r(p10)

tabstat vpk if dlne>de75&dlne~=., stat(var mean)
tabstat vpk if dlne<=de25&dlne~=., stat(var mean)


