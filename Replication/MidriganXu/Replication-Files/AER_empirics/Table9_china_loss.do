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
**generate De
gen De=tfp-L.tfp

*take the mean of tfp as z
egen zp=mean(tfp), by (plant)

egen nobs=count(year), by (plant)
gen e=tfp-zp

**identify rho, eps with fixed effect regression
xtreg tfp L.tfp, fe
**FDIV
ivregress 2sls d.tfp (Ld.tfp=LL.tfp), cluster(plant)
predict eps, res
sum eps, d
di .r(sd)/sqrt(2)

**get first best TFP for each industry
gen ee=exp(e)
egen inde=mean(ee), by (inid2 year)
gen linde=(1-0.85)*log(inde)

**************************************calculate TFP loss in the data ***********************************
gen group=cage
gen vpk=lnvp-lnk
sum vpk, d
replace vpk=. if vpk>r(p99)|vpk<r(p1)

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

**residual of group effect
gen vpkr=vpk-avpk
gen distp1r=exp(e)*(exp(vpkr)^(-(1-mal2)*0.85/(1-0.85)))
gen distp2r=exp(e)*(exp(vpkr)^((mal2*0.85-1)/(1-0.85)))

egen ydistp1=mean(distp1r), by (inid2 year group)
egen ydistp2=mean(distp2r), by (inid2 year group)

egen firstp=mean(exp(e)), by (inid2 year group)
gen yfirstp=(1-0.85)*ln(firstp)

gen gap=(1-mal2*0.85)*ln(ydistp1)-(1-mal2)*0.85*ln(ydistp2)-yfirstp

preserve
collapse (mean) gap (sum) v, by (inid2 year group)
table group year [aw=v], c(mean gap)
restore


gen diste=exp(e*(1-0.85)/(1-0.85*mal2))
egen ydiste=mean(diste), by (year inid2)
egen tfirste=mean(exp(e)), by (year inid2)
gen gape=(1-mal2*0.85)*ln(ydiste)-(1-0.85)*ln(tfirste)

*table inid2 year, c(mean gape)
preserve
collapse (mean) gape (sum) v, by (inid2 year)
tabstat gape [aw=v], by (year)
restore
