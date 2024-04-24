
clear all
set mem 200m
set more off

clear all
set mem 200m
set more off

* Change to corresponding folder

use appendix-dataset.dta

gen trade_share  = trade_goods/absorption_dgen MP_share     = sales_MandA/gross_prod_nonfin_d

gen log_share_MP = log(MP_share)     
gen log_trade_sh = log(trade_share)

gen log_M       = log(num_aff_MandA)
gen log_avg_MP  = log(MP_share)  - log(num_aff_MandA)

gen log_sales  = log(sales_MandA)
gen log_stocks = log(stocks)

gen log_dist    = log(distance)
gen log_gdp_d   = log(gdp_d)
gen log_Y_d     = log(gross_prod_nonfin_d)

gen d_fdi = 1 if sales_MandA>0
replace d_fdi = 0 if sales_MandA==0

gen d_trade = 1 if trade_goods>0
replace d_trade = 0 if trade_goods==0


tab ISO_d, gen(DD)
tab ISO_o, gen(OO)


*** GRAVITY REGRESSIONS*******************

reg d_fdi        log_dist cborder colony clanguage DD* OO*, noconstant robust
reg d_trade      log_dist cborder colony clanguage DD* OO*, noconstant robust

reg log_trade_sh log_dist cborder colony clanguage DD* OO*, noconstant robust
reg log_share_MP log_dist cborder colony clanguage DD* OO*, noconstant robust

reg log_M        log_dist cborder colony clanguage DD* OO*, noconstant robust
reg log_avg_MP   log_dist cborder colony clanguage DD* OO*, noconstant robust



***Stats for Figure 1

gen log_ekk_var = log_M - log(MP_share)

reg log_M       log_Y_d,constant robust 

reg log_ekk_var log_Y_d    ,constant robust
reg log_ekk_var log_Y_d OO*,noconstant robust


*** EKK-type regression
reg log_M  log_sales DD* OO*, noconstant robust


*****MP on FDI STOCKS (reported in "Data")
reg log_sales log_stocks DD* OO*,constant robust


* Generate country-level domestic shares
egen inward_MP_share = sum(MP_share),by(ISO_d)
egen import_share    = sum(trade_share),by(ISO_d)gen dom_share_MP    = 1 - inward_MP_share
gen dom_share_trade = 1 - import_share

egen imports   = sum(trade_goods), by(ISO_d)
egen inward_mp = sum(sales_MandA), by(ISO_d)

collapse (mean) imports inward_mp  absorption_d dom_share_MP dom_share_trade gross_prod_nonfin_d rgdpl_d gdp_d, by(ISO_d)

gen log_dom_share_trade = log(dom_share_trade)
gen log_dom_share_mp    = log(dom_share_MP) 
gen log_gross_prod 		= log(gross_prod_nonfin_d)
gen log_rgdpl      		= log(rgdpl_d)
gen log_absorption 		= log(absorption_d)


***** Regressions for the Gains from trade and MP

reg log_dom_share_trade log_rgdpl log_absorption, robust
reg log_dom_share_trade log_rgdpl log_gross_prod, robust


****Construct country-level data**********

keep imports inward_mp gdp_d ISO_d

rename gdp_d gdp
rename ISO_d ISO
save current-account-variables.dta,replace
clear

use appendix-dataset.dta
egen exports    = sum(trade_goods), by(ISO_o)
egen outward_mp = sum(sales_MandA), by(ISO_o)

collapse (mean) exports outward_mp, by(ISO_o)
rename ISO_o ISO
sort ISO
save temp.dta,replace
use current-account-variables.dta
sort ISO
merge ISO using temp.dta
drop _merge

* Profits share
gen eta = 0.1875

* Current account 
gen ca_trade = (exports - imports)/gdp
gen ca_mp    = eta*(outward_mp - inward_mp)/gdp
gen ca       = ca_mp + ca_trade


save current-account-variables.dta,replace
clear


