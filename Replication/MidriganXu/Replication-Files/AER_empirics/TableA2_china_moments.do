

clear

clear matrix

set mem 2g

cd "C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\china_data"

use china_data_clean

*capture log close

*capture log using china_inv_moments_final.smcl, replace

**0.  Number of Observations, and Unique Plants
rename id plant
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

gen cage=1 if age>=0&age<=5
replace cage=2 if age>5&age<=10
replace cage=3 if age>10&age~=.

*0.2 Define key variables
gen vm=csales-wage
replace vm=. if vm<=0

gen am=vm/output
replace am=. if am>1|am<0
egen mam=mean(am), by (inid2)

gen v=(output-vm)/cpi
gen lnv=ln(v)

gen qk=fixasset if fixasset>=0
gen ql=(wage+welfare)/cpi if wkr>=0
gen lnk=ln(qk)
gen lnl=ln(ql)


gen lnq=ln(output/cpi)
gen lnm=ln(vm/cpi)

**1.  Size Distribution and Concentration

*1.1 Dispersion of Lnv, Lnk, Lnl
sum lnv, d
matrix mdispv=[r(sd)]

sum lnk, d
matrix mdispk=[r(sd)]

sum lnl, d
matrix mdispl=[r(sd)]

tabstat lnv, by (cage) stat (var) save
matrix mdispva=[r(Stat1),r(Stat2),r(Stat3)]

tabstat lnl, by (cage) stat (var) save
matrix mdispla=[r(Stat1),r(Stat2),r(Stat3)]

tabstat lnk, by (cage) stat (var) save
matrix mdispka=[r(Stat1),r(Stat2),r(Stat3)]

*1.2 Concentration
egen yobs=count(year), by (year)
egen yaobs=count(year), by (year cage)
replace yaobs=yaobs/yobs
tabstat yaobs, by (cage) stat (mean) save
matrix ashare=[r(Stat1),r(Stat2),r(Stat3)]

gsort year -v
by year: gen rankv=[_n] 
egen tv=sum(v), by (year)

gen rankv1=1 if rankv/yobs<=0.01
gen rankv5=1 if rankv/yobs<=0.05
gen rankv10=1 if rankv/yobs<=0.1
gen rankv20=1 if rankv/yobs<=0.2

gen v0=v/tv/10
gen v1=v*rankv1/tv/10
gen v5=v*rankv5/tv/10
gen v10=v*rankv10/tv/10
gen v20=v*rankv20/tv/10
tabstat v1 v5 v10 v20, stat (sum) save
matrix vshare=r(StatTotal)
tabstat v0, by (cage) stat (sum) save
matrix vashare=[r(Stat1),r(Stat2),r(Stat3)]

gsort year -qk
by year: gen rankk=[_n] 
egen tk=sum(qk), by (year)

gen rankk1=1 if rankk/yobs<=0.01
gen rankk5=1 if rankk/yobs<=0.05
gen rankk10=1 if rankk/yobs<=0.1
gen rankk20=1 if rankk/yobs<=0.2

gen k0=qk/tk/10
gen k1=qk*rankk1/tk/10
gen k5=qk*rankk5/tk/10
gen k10=qk*rankk10/tk/10
gen k20=qk*rankk20/tk/10
tabstat k1 k5 k10 k20, stat (sum) save
matrix kshare=r(StatTotal)
tabstat k0, by (cage) stat (sum) save
matrix kashare=[r(Stat1),r(Stat2),r(Stat3)]

gsort year -ql
by year: gen rankl=[_n] 
egen tl=sum(ql), by (year)

gen rankl1=1 if rankl/yobs<=0.01
gen rankl5=1 if rankl/yobs<=0.05
gen rankl10=1 if rankl/yobs<=0.1
gen rankl20=1 if rankl/yobs<=0.2

gen l0=ql/tl/10
gen l1=ql*rankl1/tl/10
gen l5=ql*rankl5/tl/10
gen l10=ql*rankl10/tl/10
gen l20=ql*rankl20/tl/10
tabstat l1 l5 l10 l20, stat (sum) save
matrix lshare=r(StatTotal)
tabstat l0, by (cage) stat (sum) save
matrix lashare=[r(Stat1),r(Stat2),r(Stat3)]

gsort year cage -v
by year cage: gen rankav=[_n]
by year cage: egen tobsa=count(v)
by year cage: egen tva=sum(v)
gen ranka1=1 if rankav/tobsa<=0.01
gen ranka5=1 if rankav/tobsa<=0.05
gen ranka10=1 if rankav/tobsa<=0.1
gen ranka20=1 if rankav/tobsa<=0.2

gen va1=v*ranka1/tva/10
gen va5=v*ranka5/tva/10
gen va10=v*ranka10/tva/10
gen va20=v*ranka20/tva/10

tabstat va1, by (cage) stat (sum) save
matrix vashare1=[r(Stat1),r(Stat2),r(Stat3)]

gsort year cage -qk
by year cage: gen rankak=[_n]
by year cage: egen tka=sum(qk)
gen rankka1=1 if rankak/tobsa<=0.01
gen rankka5=1 if rankak/tobsa<=0.05
gen rankka10=1 if rankak/tobsa<=0.1
gen rankka20=1 if rankak/tobsa<=0.2

gen ka1=qk*rankka1/tka/10
gen ka5=qk*rankka5/tka/10
gen ka10=qk*rankka10/tka/10
gen ka20=qk*rankka20/tka/10

tabstat ka1, by (cage) stat (sum) save
matrix kashare1=[r(Stat1),r(Stat2),r(Stat3)]

gsort year cage -ql
by year cage: gen rankal=[_n]
by year cage: egen tla=sum(ql)
gen rankla1=1 if rankal/tobsa<=0.01
gen rankla5=1 if rankal/tobsa<=0.05
gen rankla10=1 if rankal/tobsa<=0.1
gen rankla20=1 if rankal/tobsa<=0.2

gen la1=ql*rankla1/tla/10
gen la5=ql*rankla5/tla/10
gen la10=ql*rankla10/tla/10
gen la20=ql*rankla20/tla/10

tabstat la1, by (cage) stat (sum) save
matrix lashare1=[r(Stat1),r(Stat2),r(Stat3)]
**2. Correlation of Y overtime
egen iplant=group(plant)

xtset iplant year

foreach var of varlist lnv lnl lnk{
gen d`var'=`var'-L.`var'
by iplant: gen `var'1=`var'[_n-1]
by iplant: gen `var'3=`var'[_n-3]
by iplant: gen `var'5=`var'[_n-5]
corr `var' `var'1
matrix corr`var'1=r(rho)
corr `var' `var'3
matrix corr`var'3=r(rho)
corr `var' `var'5
matrix corr`var'5=r(rho)
sum d`var',d
matrix sd`var'=r(sd)
}
/*gen dlnv=lnv-L.lnv
gen dlnl=lnl-L.lnl
gen dlnk=lnk-L.lnk*/
/*
by iplant: gen lnv1=lnv[_n-1]
by iplant: gen lnl1=lnl[_n-1]
by iplant: gen lnk1=lnk[_n-1]
by iplant: gen lnv3=lnv[_n-3]
by iplant: gen lnl3=lnl[_n-3]
by iplant: gen lnk3=lnk[_n-3]
by iplant: gen lnv5=lnv[_n-5]
by iplant: gen lnl5=lnl[_n-5]
by iplant: gen lnk5=lnk[_n-5]*/

/*corr lnv lnv1
matrix corrlnv1=r(rho)
corr lnv lnv3
matrix corrlnv3=r(rho)
corr lnv lnv5
matrix corrlnv5=r(rho)*/


/*sum dlnv, d
matrix sdlnv=r(sd)
sum dlnl,d
matrix sdlnl=r(sd)
sum dlnk,d
matrix sdlnk=r(sd)*/

**3. Average revenue product of K, employment to K ratio, by age
gen vpk=log(v/qk)
gen lk=log(ql/qk)
gen lnda=totdebt/totasset
replace lnda=. if lnda>1|lnda<0

egen indy=group(inid2 year)
areg vpk, absorb(indy)
predict vpkk, res
replace vpkk=vpkk+_coef[_cons]

areg lk, absorb(indy)
predict lkk, res
replace lkk=lkk+_coef[_cons]

areg lnda, absorb(indy)
predict lndap, res
replace lndap=lndap+_coef[_cons]

tabstat vpkk, by (cage) stat (mean) save
matrix mvpk=[r(Stat1),r(Stat2),r(Stat3)]
tabstat lkk, by (cage) stat (mean) save
matrix mlkk=[r(Stat1),r(Stat2),r(Stat3)]
tabstat lndap, by (cage) stat(mean) save
matrix mlndap=[r(Stat1),r(Stat2),r(Stat3)]

tabstat vpkk, by (cage) stat (var) save
matrix sdvpk=[r(Stat1), r(Stat2),r(Stat3)]
tabstat lkk, by (cage) stat (var) save
matrix sdlkk=[r(Stat1),r(Stat2),r(Stat3)]

tab cage, gen (dage)

reg vpkk dage1 dage2 lnv
matrix rvpk=[_coef[dage1],_coef[dage2]]
reg lkk dage1 dage2 lnl
matrix rlkk=[_coef[dage1],_coef[dage2]]
reg lndap dage1 dage2 lnv
matrix rlndapv=[_coef[dage1],_coef[dage2]]
reg lndap dage1 dage2 lnl
matrix rlndapl=[_coef[dage1],_coef[dage2]]

xtset iplant year
gen dvpk=vpkk-L.vpkk
gen dvpk3=vpkk-L.L.L.vpkk
gen dvpk5=vpkk-L.L.L.L.L.vpkk

gen dlkk=lkk-L.lkk

areg vpk dlnv L.vpk, absorb(indy)
matrix evpk=[_coef[dlnv],_coef[L.vpk]]
/*areg dvpk3 dlnv3, absorb(indy)
areg dvpk5 dlnv5, absorb(indy)*/

areg lkk dlnv L.lkk, absorb(indy)
matrix elkk=[_coef[dlnv],_coef[L.lkk]]


**4. Share of exits
egen txv=sum(v) if maxyear==year, by (year)
gen xshare=txv/tv
tabstat xshare, by (year) stat (mean) save
matrix mxshare=[(r(Stat1)+r(Stat2)+r(Stat3)+r(Stat4)+r(Stat5)+r(Stat6)+r(Stat7)+r(Stat8)+r(Stat9))/9]

**5. Start with Age aggregate moments - Growth rate by age

sort plant year
by plant: gen nv=v[_n+1]
by plant: gen nqk=qk[_n+1]
by plant: gen nql=ql[_n+1]

egen sv=sum(v) if nv~=., by (year cage)
egen sqk=sum(qk) if nqk~=., by (year cage)	
egen sql=sum(ql) if nql~=., by (year cage)

egen snv=sum(nv), by (year cage)
egen snqk=sum(nqk), by (year cage)
egen snql=sum(nql), by (year cage)

collapse (mean) sv snv sqk snqk sql snql (sum) v ql qk totdebt totasset, by (inid2 year cage)
gen vpk=log(v/qk)
gen lk=log(ql/qk)
gen lnda=totdebt/totasset
sum lnda, d

egen indy=group(inid2 year)
areg vpk, absorb(indy)
predict vpkk, res
replace vpkk=vpkk+_coef[_cons]

areg lk, absorb(indy)
predict lkk, res
replace lkk=lkk+_coef[_cons]

areg lnda, absorb(indy)
predict lndap, res
replace lndap=lndap+_coef[_cons]

gen dlnva=log(snv)-log(sv)
gen dlnka=log(snqk)-log(sqk)
gen dlnla=log(snql)-log(sql)

collapse (mean) vpkk lkk lndap dlnva dlnka dlnla, by (year cage)

tabstat dlnva, by (cage) stat (mean) save	
matrix adlnva=[r(Stat1)-r(Stat3),r(Stat2)-r(Stat3)]
tabstat dlnka, by (cage) stat (mean) save	
matrix adlnka=[r(Stat1)-r(Stat3),r(Stat2)-r(Stat3)]
tabstat dlnla, by (cage) stat (mean) save	
matrix adlnla=[r(Stat1)-r(Stat3),r(Stat2)-r(Stat3)]

tabstat vpkk, by (cage) stat (mean) save
matrix avpk=[r(Stat1),r(Stat2),r(Stat3)]
tabstat lkk, by (cage) stat (mean) save
matrix alkk=[r(Stat1),r(Stat2),r(Stat3)]
tabstat lndap, by (cage) stat(mean) save
matrix alndap=[r(Stat1),r(Stat2),r(Stat3)]


**summarize
matrix datam=[sdlnv'\sdlnl'\sdlnk'\mdispv'\mdispl'\mdispk'\corrlnv1'\corrlnv3'\corrlnv5'\corrlnl1'\corrlnl3'\corrlnl5'\corrlnk1'\corrlnk3'\corrlnk5'\vshare'\mxshare'\mdispva'\mdispla'\mdispka'\vashare1'\lashare1'\kashare1'\ashare'\vashare'\lashare'\kashare'\avpk'\mvpk'\alkk'\mlkk'\adlnva'\adlnla'\adlnka'\sdvpk'\sdlkk'\alndap'\mlndap'\rvpk'\rlkk'\evpk'\elkk'\rlndapv'\rlndapl']
drop year-dlnla
svmat datam
