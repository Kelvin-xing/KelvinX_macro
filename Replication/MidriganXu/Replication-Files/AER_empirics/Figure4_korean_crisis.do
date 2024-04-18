clear

clear matrix

set mem 500m

clear

clear matrix

set mem 1g

cd "C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\korean_data"

use "C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\korean_data\raw data processing\korea_data.dta"

sort busstr year
egen plant=group(busstr)

sort year inid2

keep if year<=99

sort busstr year
by busstr: egen minyear=min(year)
by busstr: egen maxyear=max(year)

gen ecap_b=(bldgby+mchby+carby)/iprice  /**book value BOY for comparison**/
gen ecap_p=(bldgby+mchby+carby)/iprice  if year==minyear 

gen einv=(bdpch+mchpch+carpch-bdsell-mcsell-casell)/iprice
gen eret=(bddep+mchdep+cardep)/iprice

by busstr: replace ecap_p=ecap_p[_n-1]+einv[_n-1]-eret[_n-1] if year~=minyear

tabstat ecap_p, by (year) stat (min p1 p5 p25 p50 p75 p95 p99 max)
tabstat ecap_b, by (year) stat (min p1 p5 p25 p50 p75 p95 p99 max)

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
egen mam=mean(am), by (inid2)

*labor input 
gen ql=(twage+wfare)/cpi

*capital and investment
gen qk=ecap_p +rent/cpi/uk
gen qi=einv

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
table year flagdrop, c(sum v)

*********************************1. report aggregate statistics with/without outliers**********************************
drop if dumdrop==1   /**negative values, but at plant-year level,  other obs of plant remained**/
tab year
tabstat v ql qk qe, by (year) stat (sum)

drop if flagdrop==1
tab year
tabstat v ql qk qe, by (year) stat (sum)

unique plant


*****FOR MICRO MOMENTS, GET RID OF PLANTS ENTERING AT 99 (MIS-REPORTED CAPITAL MEASURE)

*********************************2. define age variable and construct plant-level  tfp_it*******************************
/*define age/cohort*/
replace yborn=. if yborn<0|yborn>99

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

gen tfp1=lnv-mal*lnl-0.85*(1-mal)*lnk
gen tfpp1=lnvp-mal*lnl-0.85*(1-mal)*lnk

gen mal2=mal*(0.66/0.52)
gen tfp=(lnvp-mal2*0.85*lnl-(1-mal2)*0.85*lnk)/(1-0.85)
sum tfp, d
replace tfp=. if tfp>r(p99)|tfp<r(p1)

xtset plant year
**generate De
gen De=tfp-L.tfp

gen dcrisis=1 if year>=97
replace dcrisis=0 if dcrisis==.
gen ctfp=dcrisis*L.tfp


*separate tranistory shock vs permanent one

*take the mean of tfp as z
egen zp=mean(tfp), by (plant)

egen nobs=count(year), by (plant)
gen e=tfp-zp 

**get first best TFP for each industry
gen ee=exp(e)
egen inde=mean(ee), by (inid2 year)
gen linde=(1-0.85)*log(inde)


***define which group we want to use
gen group=cage

****************************************simple moments based on groups *****************
gen vpk=lnvp-lnk
sum vpk, d
replace vpk=. if vpk>r(p99)|vpk<r(p1)

xtset plant year
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

foreach y in 92 93 94 95 96 97 98 99{
areg vpk e L.vpk if year==`y'&cage==1, absorb(inid2)
areg vpk e L.vpk if year==`y'&cage==3, absorb(inid2)
}
*************************************elast. D(Y or K) on De  ******************
tab group, gen (dgroup)
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

matrix datam=[mvpk'\mvpka'\svpk'\svpka'\sdlnv'\sdlnva'\sdlnk'\sdlnka'\sdtfp'\sdtfpa'\redlnv'\redlnva'\redlnk'\redlnka'\redvpk'\redvpka'\rvdlnk'\rvdlnka'\rvdvpk'\rvdvpka']

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
collapse (mean) tgap linde (count) nplant=e (sum) v, by (inid2 year)

gen nplantx=ln(nplant)*(1-0.85)

table inid2 year, c(mean linde)
table inid2 year, c(mean nplantx)
table inid2 year, c(mean nplant)
table inid2 year, c(mean tgap)

table year [aw=v], c(mean tgap)
table year [aw=v], c(mean linde)
table year [aw=v], c(mean nplantx)
table year [aw=v], c(mean nplant)
restore

/**due to group**/
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

**loss by age groups
gen vpkr=vpk-avpk
gen distp1r=exp(e)*(exp(vpkr)^(-(1-mal2)*0.85/(1-0.85)))
gen distp2r=exp(e)*(exp(vpkr)^((mal2*0.85-1)/(1-0.85)))

egen ydistp1=mean(distp1r), by (inid2 year group)
egen ydistp2=mean(distp2r), by (inid2 year group)

egen firstp=mean(exp(e)), by (inid2 year group)
gen yfirstp=(1-0.85)*ln(firstp)

gen gap=(1-mal2*0.85)*ln(ydistp1)-(1-mal2)*0.85*ln(ydistp2)-yfirstp

preserve
collapse (mean) gap linde (count) nplant=e (sum) v, by (inid2 year group)

replace nplant=ln(nplant)*(1-0.85)
table group year [aw=v], c(mean linde)
table group year [aw=v], c(mean nplant)
table group year [aw=v], c(mean gap)
restore


gen diste=exp(e*(1-0.85)/(1-0.85*mal2))
egen ydiste=mean(diste), by (year inid2)
egen tfirste=mean(exp(e)), by (year inid2)
gen gape=(1-mal2*0.85)*ln(ydiste)-(1-0.85)*ln(tfirste)

preserve
collapse (mean) gape (sum) v, by (inid2 year)
tabstat gape [aw=v], by (year)
restore

