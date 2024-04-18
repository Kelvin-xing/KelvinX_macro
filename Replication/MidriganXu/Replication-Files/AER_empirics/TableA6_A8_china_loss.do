
/***We illustrate using China, Korea and Colombia are constructed in a very similar way***/

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

*0.2 Start to define alternative group of plants

**based on age
gen cage=1 if age>=0&age<=5
replace cage=2 if age>5&age<=10
replace cage=3 if age>10&age~=.

**based on industry DEF 
gen cdef=1 if inid2==13|inid2==14|inid2==15|inid2==16|inid2==18|inid2==19|inid2==22|inid2==23|inid2==24|inid2==31|inid2==32|inid2==33
replace cdef=2 if inid2==17|inid2==20|inid2==21|inid2==25|inid2==26|inid2==28|inid2==29|inid2==34|inid2==37
replace cdef=3 if inid2==27|inid2==30|inid2==31|inid2==35|inid2==36|inid2==39|inid2==40|inid2==41

**based on debt-to-capital ratio
replace fixasset=. if fixasset<0
gen ndebt=totdebt-(totasset-fixasset)
gen dk=ndebt/fixasset
gen cdebt=1 if dk<=0&dk~=.
replace cdebt=2 if dk>0&dk<.4&dk~=.
replace cdebt=3 if dk>=.4&dk~=.

**based on ownership
replace owner=3 if owner==4

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

xtset plant year

*take the mean of tfp as z
egen zp=mean(tfp), by (plant)

egen nobs=count(year), by (plant)
gen e=tfp-zp

**get first best TFP for each industry
gen ee=exp(e)
egen inde=mean(ee), by (inid2 year)
gen linde=(1-0.85)*log(inde)


***define which group we want to use
*gen group=cdef
*gen group=cdebt
gen group=owner
tab group, gen (dgroup)

****************************************simple moments based on groups *****************
gen vpk=lnvp-lnk
sum vpk, d
replace vpk=. if vpk>r(p99)|vpk<r(p1)

xtset plant year

reg vpk dgroup2 dgroup3 i.year 

areg vpk, absorb(indy)
predict shifter, d
gen vpkk=vpk-shifter

robvar vpkk, by (group) 

**************************************calculate TFP loss in the data ***********************************

**total misallocation due to capital
gen distp1=exp(e)*(exp(vpk)^(-(1-mal2)*0.85/(1-0.85)))
gen distp2=exp(e)*(exp(vpk)^((mal2*0.85-1)/(1-0.85)))

egen tydistp1=mean(distp1), by (year inid2)
egen tydistp2=mean(distp2), by (year inid2)
egen tfirstp=mean(exp(e)), by (year inid2)
gen tyfirstp=(1-0.85)*ln(tfirstp)
gen tgap=(1-mal2*0.85)*ln(tydistp1)-(1-mal2)*0.85*ln(tydistp2)-tyfirstp

preserve
collapse (mean) tgap (sum) v, by (inid2 year)
table year [aw=v], c(mean tgap)
restore

**due to group
egen avpk=mean(vpk), by (year inid2 group)
gen distp1a=exp(e)*(exp(avpk)^(-(1-mal2)*0.85/(1-0.85)))
gen distp2a=exp(e)*(exp(avpk)^((mal2*0.85-1)/(1-0.85)))

egen ydistp1a=mean(distp1a), by (year inid2 )
egen ydistp2a=mean(distp2a), by (year inid2)
gen gapa=(1-mal2*0.85)*ln(ydistp1a)-(1-mal2)*0.85*ln(ydistp2a)-tyfirstp

preserve
collapse (mean) gapa (sum) v, by (inid2 year)
tabstat gapa [aw=v], by (year)
restore

