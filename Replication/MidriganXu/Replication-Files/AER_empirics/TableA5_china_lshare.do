clear

clear matrix

set mem 1g

cd "C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\china_data"

use china_data_clean

**0.  Number of Observations, and Unique Plants
egen plant=group(id)
gen inid2=floor(cic_adj/100)

egen nyear=count(year), by (plant)
drop if nyear>10
sort plant year
by plant: gen error_y=1 if year==year[_n-1]
drop if error_y==1

tab year
unique plant

egen minyear=min(year), by (plant)
egen maxyear=max(year), by (plant)

*0.1 Define age/cohort
gen age=year-cohort
**check miss reporting
replace age=. if age<0
**
replace age=. if age==0&cohort>minyear
drop if age==.

*0.3 Define key variables
gen vm=csales-wage
replace vm=. if vm<=0

gen am=vm/output
replace am=. if am>1|am<0
egen mam=mean(am), by (inid2)

gen v=(output-vm)/cpi
replace v=. if v<=0
gen lnv=log(v)

gen qk=fixasset if fixasset>=0
gen ql=(wage+welfare)/cpi if wkr>=0
gen lnk=log(qk)
gen lnl=log(ql)

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

gen mal2=mal*(0.66/0.32)
gen tfp2=lnv-mal2*0.85*lnl-(1-mal2)*0.85*lnk
gen tfp=(lnvp-mal2*0.85*lnl-(1-mal2)*0.85*lnk)/(1-0.85)

sum tfp, d
replace tfp=. if tfp>r(p99)|tfp<r(p1)

tabstat mal2, by (inid2)
