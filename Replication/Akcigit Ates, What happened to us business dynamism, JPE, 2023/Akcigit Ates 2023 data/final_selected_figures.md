---
title: "Final Selected Figures"
author: "Zhiyu Fu"
date: "1/27/2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


This file are designed to be the most bottom-up and self-contained file that can reproduce all the results in the paper from raw data.

Links to raw data:

[patent.tsv](http://www.patentsview.org/data/20171226/patent.tsv.zip)
[patent_assignee.tsv](http://s3.amazonaws.com/data.patentsview.org/20181127/download/patent_assignee.tsv.zip)
[application.tsv](http://s3.amazonaws.com/data.patentsview.org/20181127/download/application.tsv.zip)
[assignee.tsv](http://s3.amazonaws.com/data.patentsview.org/20181127/download/assignee.tsv.zip)
[apat63_99.csv](https://www-nber-org.proxy.uchicago.edu/patents/apat63_99.zip)
[reassignment datasets](https://bulkdata.uspto.gov/data/patent/assignment/economics/2017/csv.zip)
[uspatentcitation.tsv](http://www.patentsview.org/data/20171226/uspatentcitation.tsv.zip) 
[MaintFeeEvents](https://bulkdata.uspto.gov/data/patent/maintenancefee/MaintFeeEvents.zip)

Claim length: TBS.

```{r}
library(tidyverse)
library(data.table)
library(knitr)
library(foreign)
```


# Average filing patents per firm

Links to backbone datasets:



```{r}
patent_assignee = read_delim("../raw_data/patent_assignee.tsv", delim = '\t',
                             col_types = cols_only(patent_id = 'c', assignee_id = 'c')) %>% data.table()
application = read_delim("../raw_data/application.tsv", delim = '\t',
          col_types = cols_only(id = 'c', patent_id = 'c', date = 'D', number = 'c')) %>%
  rename(application_id = id, application_number = number, appdate = date) %>% data.table()
patent = read_delim("../raw_data/patent.tsv", delim = '\t',
           col_types = cols_only(id = 'c', type = 'c', date = 'D')) %>%
  rename(patent_id = id, gdate = date, patent_type = type) %>% data.table()
assignee = read_delim("../raw_data/assignee.tsv", delim = '\t', col_types = cols_only(
    id = 'c',
    type = 'i',
  organization = 'c'
)) %>% rename(assignee_id = id, assignee_type = type) %>% data.table()
full_patent = right_join(patent_assignee, patent)
full_patent = left_join(full_patent, application)
full_patent = left_join(full_patent, assignee)
save(full_patent, file = "full_patent.RData")
```

```{r}
full_patent = data.table(full_patent)
usfirm_patent = full_patent[patent_type == "utility" & assignee_type == 2] %>%
  select(-patent_type, -application_id, assignee_type) %>%
  mutate(appyear = year(appdate), gyear = year(gdate)) %>% data.table()
usfirm_patent = usfirm_patent[order(patent_id, assignee_id)]
usfirm_patent = usfirm_patent[!duplicated(patent_id)]
```

```{r}
# Calculate size for each active firm across time
pat_assignee = usfirm_patent[!is.na(appyear),.(assignee_id, appyear)]
pat_assignee_crossed = inner_join(pat_assignee[,.(assignee_id, year = appyear)],
                          unique(pat_assignee)) %>% data.table()
assignee_year_size = pat_assignee_crossed[,.(
  cumN = sum(year<=appyear),
  N5yr = sum(year<=appyear & (year > (appyear-5))),
  N1yr = sum(year==appyear)
),by = .(assignee_id, appyear)]
usfirm_patent = left_join(usfirm_patent, assignee_year_size, by = c('assignee_id', 'appyear')) %>%
  data.table()
remove(pat_assignee_crossed, pat_assignee)
usfirm_patent = usfirm_patent[order(appyear, assignee_id)]
save(usfirm_patent, file = 'usfirm_patent.RData')
save(assignee_year_size, file = "assignee_year_size.RData")
# load('usfirm_patent.RData')
```

```{r}
# Calculate the firm size distribution (measured by patents) across time
pat_assignee = usfirm_patent[!is.na(appyear),.(assignee_id, appyear)]
series_perc = list()
pat_distr_list = list()
for (year in 1975:2017)
{
  pat_assignee[,flag1yr:= appyear == year]
  pat_assignee[, flag5yr:= appyear <= year & appyear > year - 5]
  distr1yr = pat_assignee[,.(npat1yr = sum(flag1yr)), by = assignee_id] %>% 
    unique() %>% filter(npat1yr>=1) 
  distr5yr = pat_assignee[,.(npat5yr = sum(flag5yr)), by = assignee_id] %>% 
    unique() %>% filter(npat5yr>=1) 
  # whether use npat1yr>=1 or npat5yr>=1 depends on situations. If our target firm may not innovate in the targeted year, we should filter by npat5yr. use distr5yr2 to denote sample filtered by npat1yr.
  distr5yr2 = pat_assignee[,.(npat5yr = sum(flag5yr), npat1yr = sum(flag1yr)), by = assignee_id] %>% 
    unique() %>% filter(npat1yr>=1) 
  perc5yr = distr5yr %>% 
    summarize(
              q5_5yr = quantile(npat5yr, 0.05),
              q50_5yr = quantile(npat5yr, 0.50),
              q90_5yr = quantile(npat5yr, 0.90),
              q99_5yr = quantile(npat5yr, 0.99),
              mean5yr = mean(npat5yr),
              meanlog5yr = mean(log(npat5yr)),
              sd5yr = sd(npat5yr),
              n_firm5yr = length(assignee_id),
              year = year)
  perc5yr2 = distr5yr2 %>% 
    summarize(
              q5_5yr2 = quantile(npat5yr, 0.05),
              q50_5yr2 = quantile(npat5yr, 0.50),
              q90_5yr2 = quantile(npat5yr, 0.90),
              q99_5yr2 = quantile(npat5yr, 0.99),
              mean5yr2 = mean(npat5yr),
              meanlog5yr2 = mean(log(npat5yr)),
              sd5yr2 = sd(npat5yr),
              n_firm5yr2 = length(assignee_id),
              year = year)
  perc1yr = distr1yr %>%
    summarize(
              q5_1yr = quantile(npat1yr, 0.05),
              q50_1yr = quantile(npat1yr, 0.50),
              q90_1yr = quantile(npat1yr, 0.90),
              q99_1yr = quantile(npat1yr, 0.99),
              mean1yr = mean(npat1yr), 
              meanlog1yr = mean(log(npat1yr)),
              sd1yr = sd(npat1yr),
              n_firm1yr = length(assignee_id),
              year = year
              )
  distr = list(distr1yr, distr5yr, distr5yr2)
  names(distr) = c('1yr', '5yr', '5yr2')
  perc = inner_join(perc1yr, perc5yr, by= 'year') %>% inner_join(perc5yr2, by = 'year')
  series_perc = append(series_perc, list(perc))
  pat_distr_list = append(pat_distr_list, list(distr))
}
names(pat_distr_list) = 1975:2016
pat_distr_perfirm = rbindlist(series_perc) %>% data.table()
remove(perc1yr, perc5yr, perc, distr1yr, distr5yr, distr, distr5yr2, perc5yr2)
save(pat_distr_perfirm, file = "pat_distr_perfirm.RData")
# load('pat_distr_perfirm.RData')
```



```{r}
nber_pat = read_csv("../raw_data/apat63_99.csv", 
                    col_types = cols_only(ASSIGNEE = 'c', PATENT = 'c', APPYEAR = 'i', 
                                          GYEAR = 'c', ASSCODE = 'c')) %>% data.table()
colnames(nber_pat) = sapply(colnames(nber_pat), str_to_lower)
nber_uspto_matched = nber_pat[!is.na(assignee)] %>% select(nber_assignee = assignee, patent_id = patent, 
                                                  nber_appyear = appyear, 
                                                  nber_gyear = gyear, 
                                                  nber_asscode = asscode) %>%
  mutate(patent_id = as.character(patent_id)) %>%
  inner_join(usfirm_patent) %>% data.table()

# Check the consistency of two dataset.
nber_uspto_matched[, sum(is.na(nber_appyear) |  is.na(appyear) |nber_appyear!=appyear)]
nber_uspto_matched[, sum(is.na(nber_gyear) |  is.na(gyear) |nber_gyear!=gyear)]
nber_uspto_matched[,table(nber_asscode)]

# quite consistent
crosswalk = nber_uspto_matched[, .(nber_assignee, assignee_id)] %>% unique()
crosswalk[,Nnberassg := .N, by = nber_assignee]
crosswalk[,Nusptoassg := .N, by = assignee_id]

# Notice the mapping is complicated.
crosswalk[,mean(Nnberassg>1)]

# We want to map nber assignee to uspto assignee. Thus, let's drop duplicated nber assignee.
crosswalk = crosswalk[!duplicated(nber_assignee)]
nber_pat = merge(crosswalk, nber_pat, by.x = "nber_assignee", by.y = "assignee", all.y = T)
usfirm_nber_pat_pre1975 = nber_pat[asscode==2][gyear<1975]
usfirm_nber_pat_pre1975[, mean(is.na(assignee_id))]
usfirm_nber_pat_pre1975[, mean(is.na(nber_assignee))]
usfirm_nber_pat_pre1975 = usfirm_nber_pat_pre1975[,.(patent_id = patent, assignee_id, gyear, appyear, assignee_type = asscode)] %>%
  mutate(patent_id = as.character(patent_id)) 

usfirm_uspto_patent_post1975 = usfirm_patent[gyear>=1975]
usfirm_patent_1963_2017 = rbind(usfirm_uspto_patent_post1975, usfirm_nber_pat_pre1975, fill = T)
save(usfirm_patent_1963_2017, file = "usfirm_patent_1963_2017.RData")
```


```{r}
usfirm_patent_1963_2017[,appyear := as.integer(appyear)]
usfirm_patent_1963_2017 = usfirm_patent_1963_2017[order(assignee_id, appyear)]
usfirm_patent_1963_2017[,by_entrant:=!duplicated(assignee_id)]
usfirm_patent_1963_2017[,.(n_pat_by_entrant = sum(by_entrant), n_pat = .N), by = appyear] %>%  
  filter(appyear %in% 1975:2014) %>% ggplot(aes(appyear, n_pat_by_entrant/n_pat)) + 
  geom_line() + geom_point()+
  scale_y_continuous(name = '', labels = scales::percent) +
  scale_x_continuous(name = "Year of Application") 
ggsave("entrant_share.png", width = 10, height = 6)
usfirm_patent_1963_2017[,.(entrant_share = sum(by_entrant)/.N), by = appyear] %>%  
  filter(appyear %in% 1975:2014) %>% write.dta(., "dta_for_figures/entrant_share.dta")
```



```{r}
pat_distr_perfirm = pat_distr_perfirm %>% data.table()
usfirm_patent = pat_distr_perfirm[,.(q99_5yr2, appyear = year)] %>% inner_join(usfirm_patent) %>%
  data.table()
usfirm_patent[,top1:=N5yr>=q99_5yr2]
# usfirm_patent[is.na(top1)] # no missing
usfirm_patent[appyear <=2014,mean(top1),by = appyear] %>% ggplot(aes(appyear, V1)) + geom_line() + geom_point()+labs(title = "Top 1% Patenting Share") +
  scale_y_continuous(name = '', labels = scales::percent) +
  scale_x_continuous(name = "Year of Application") 

ggsave("top1_patent_share.png", width = 10, height = 6)
usfirm_patent[appyear <=2014,.(top1_patent_share = mean(top1)),by = appyear] %>%
  write.dta("dta_for_figures/top1_patent_share.dta")
```


```{r}
usfirm_patent[appyear <=2014,.(patN = .N), by = appyear] %>% inner_join(pat_distr_perfirm, by = c('appyear' = 'year')) %>% 
  ggplot(aes(appyear, patN/n_firm1yr)) + geom_line() + geom_point() + labs(title = "Patents per Innovating Firms")
ggsave("patents_per_firm.png", width = 10, height = 6)
usfirm_patent[appyear <=2014,.(patN = .N), by = appyear] %>% inner_join(pat_distr_perfirm, by = c('appyear' = 'year')) %>%
  mutate(mean_size = patN/n_firm1yr) %>% select(appyear, mean_size) %>% write.dta(., "dta_for_figures/patents_per_firm.dta")
```

# Reassignment


```{r}
rm(list=ls())

assg_convey = read_csv("../raw_data/reassignment/assignment_conveyance.csv", col_types = cols_only(rf_id = 'i',
                                         employer_assign = 'i',
                                         convey_ty = 'c'))
assg = read_csv("../raw_data/reassignment/assignment.csv",
                   col_types = cols_only(rf_id = 'i',
                                         record_dt = 'D',
                                         convey_text = 'c'))
assignee = read_csv("../raw_data/reassignment/assignee.csv",
                    col_types = cols_only(rf_id = 'i',
                                          ee_name = 'c'))
assignor = read_csv("../raw_data/reassignment/assignor.csv",
                   col_types = cols_only(rf_id = 'i',
                                         exec_dt = 'D',
                                         or_name = 'c'))
doc = read_csv("../raw_data/reassignment/documentid.csv",
               col_types = cols_only(rf_id = 'i',
                                     grant_doc_num = 'c'))
doc_admin = read_csv("../raw_data/reassignment/documentid_admin.csv",
                     col_types = cols_only(rf_id = 'i',
                                           grant_doc_num = 'c',
                                           error = 'c'))
assg_convey = assg_convey %>% data.table()
assg = data.table(assg)
assignee = assignee %>% data.table()
assignor = assignor %>% data.table()
doc = doc %>% data.table()
doc_admin = data.table(doc_admin)
# save(assg, assg_convey, assignee, assignor, doc_admin, file = "../raw_data/reassign.RData")
```

```{r}

load(file = "../raw_data/reassign.RData")
load("usfirm_patent.RData")
load("pat_distr_perfirm.RData")
load("assignee_year_size.RData")



usfirm_patent = usfirm_patent[order(assignee_id, appdate)]
usfirm_patent[,entrant_patent := !duplicated(assignee_id)]
usfirm_patent = usfirm_patent[order(assignee_id, appdate)]
usfirm_patent = usfirm_patent[,enter_year:=max(entrant_patent), by = .(assignee_id, appyear)]

# Identify those patents filed by top 1%/10% firms (unused for the current version)
usfirm_patent = pat_distr_perfirm[,.(q99_5yr2, q90_5yr2, appyear = year)] %>% 
  left_join(usfirm_patent, .) %>% data.table()
usfirm_patent[,filedby10perc2:=N5yr>q90_5yr2]
usfirm_patent[,filedby1perc2:=N5yr>q99_5yr2]



# Identify the top firms at the end of 2017 (unused for the current version)
assignee_year_size = assignee_year_size[order(-cumN, assignee_id)]
assignee_cumN = assignee_year_size[!duplicated(assignee_id)]
assignee_cumN[,top01allyear:= ecdf(cumN)(cumN)>0.999]
assignee_cumN[,top1allyear:= ecdf(cumN)(cumN)>0.99]
assignee_cumN[,top10allyear:= ecdf(cumN)(cumN)>0.9]
```



```{r}
# Diambiguition algorithm

require(data.table)
require(tidyverse)

keywords = c("AB","AG","BV","CENTER","CO","COMPANY",
             "COMPANIES","CORP","CORPORATION","DIV",
             "GMBH","GROUP","INC","INCORPORATED","KG",
             "LC","LIMITED","LIMITED PARTNERSHIP","LLC","LP","LTD","NV",
             "PLC","SA","SARL","SNC","SPA","SRL","TRUST","USA","KABUSHIKI",
             "KAISHA","AKTIENGESELLSCHAFT","AKTIEBOLAG","SE","CORPORATIN",
             "CORPORATON","TRUST","GROUP","GRP","HLDGS","HOLDINGS","COMM",
             "INDS","HLDG","TECH","GAISHA")
keywords_plus = c("AMERICAN","INDUSTRY","INDUSTRIES","INTERNATIONAL", "BUSINESS", "ADMINISTRATION", "EQUIPMENT", "PRODUCTS")
all_keywords = c(keywords, keywords_plus)
firm_identify = function(name, keywords)
{
  identifier1 = str_c(sapply(keywords, function (x) str_c(' ',x, ' ')))
  identifier2 = str_c(sapply(keywords, function (x) str_c('^',x, ' ')))
  identifier3 = str_c(sapply(keywords, function (x) str_c(' ',x, '$')))
  identifier = str_c(c(identifier1, identifier2, identifier3), collapse = '|')
  newname = str_to_upper(name)
  newname = str_replace_all(newname, '[^A-Z ]', ' ')
  str_detect(newname, identifier)
}

name_matching = function(name, keywords){
  identifier1 = str_c(sapply(keywords, function (x) str_c(' ',x, ' ')))
  identifier2 = str_c(sapply(keywords, function (x) str_c('^',x, ' ')))
  identifier3 = str_c(sapply(keywords, function (x) str_c(' ',x, '$')))
  identifier = str_c(c(identifier1, identifier2, identifier3), collapse = '|')
  name = str_to_upper(name)
  newname = str_replace(name, ',.*', '')
  newname = str_replace(newname, ';.*', '')
  newname = str_replace(newname, '^THE ', '')
  newname = str_replace_all(newname, '[^A-Z ]', ' ')
  newname = str_replace_all(newname, identifier, ' ')
  newname = str_replace_all(newname, identifier, ' ')
  newname = str_replace_all(newname, ' ', '')
  protect = str_length(newname) == 0
  newname[protect] = name[protect]
  newname[protect] = str_replace(newname[protect], '^THE ', '')
  newname[protect] = str_replace_all(newname[protect], '[^A-Z ]', ' ')
  newname[protect] = str_replace_all(newname[protect], identifier, ' ')
  newname[protect] = str_replace_all(newname[protect], identifier, ' ')
  newname = str_replace_all(newname, ' ', '')
  newname[newname == ''] = name[newname == '']
  newname
}
```

```{r}
assignee_name = usfirm_patent[!is.na(organization),.(assignee_id, organization)] %>% unique()
assignee_name[,org_dismbg:=name_matching(organization, keywords)]
usfirm_patent = assignee_name[,.(assignee_id, org_dismbg)][usfirm_patent, on = "assignee_id"]
```


```{r}
# Merge reassignment datasets
doc_admin = doc_admin[error == 'none']
assg_convey = assg_convey %>% data.table()
setkey(assg, 'rf_id')
setkey(assg_convey, 'rf_id')
setkey(assignee, 'rf_id')
setkey(assignor, 'rf_id')
setkey(doc_admin, 'rf_id')
assg = assg[assg_convey]
assg[,record_year:=year(record_dt)]
assg_doc = assg[doc_admin[,.(rf_id, grant_doc_num)], nomatch = 0]
assg_doc = assg_doc[!is.na(grant_doc_num)] %>% rename(patent_id = grant_doc_num)
```

```{r}
# disambiguate names of assignees and assignors.

eename = assignee[!is.na(ee_name),.(ee_name)] %>% unique()
eename[,ee_firm:=firm_identify(ee_name, all_keywords)]
eename[,ee_dismbg:=name_matching(ee_name, keywords)]

orname = assignor[!is.na(or_name),.(or_name)] %>% unique()
orname[,or_firm:=firm_identify(or_name, all_keywords)]
orname[,or_dismbg:=name_matching(or_name, keywords)]

save(orname, eename, file = 'or_ee_name_dismbg.RData')
# load('or_ee_name_dismbg.RData')
assg_or_ee = assignor[assignee[assg, on = 'rf_id'], on = 'rf_id']
assg_or_ee = orname[eename[assg_or_ee, on = 'ee_name'], on = 'or_name']
remove(assignee, assignor)
```


Notice that there are a considerable proportion of "employer assignment" (Marco, 2015) in the dataset. Since we are only interested in inter-firm reassignment, we can filter them out by only keeping the transactions where assignors are organizations, identified by our algorithm. We did this as follows:

```{r}
assg_or_ee[,entry_N := .N, by = rf_id]
assg_or_ee[,rfid_or_firm := or_firm, by = rf_id]
assg_or_ee[is.na(rfid_or_firm), rfid_or_firm:=F]
assg_or_ee[,rfid_or_firm:= as.logical(max(rfid_or_firm)), by = rf_id]
# as long as there is at least one assignor is an organization, I'll keep it in the sample.
# most of assignees are firms, and those marked as nonfirms are mostly due to the inaccuracy of the algorithm.
assg_or_ee[is.na(exec_dt), exec_dt:=record_dt]
assg_or_ee[,exec_year:=year(exec_dt)]
assg_or_ee[,exec_year:=max(exec_year), by = rf_id] # use the latest year as the exec_year.
assg_or_ee[,mycoo :=  rfid_or_firm==T & (convey_ty == 'assignment' | convey_ty == 'merger')]
assg_coo = assg_or_ee[,.(rfid_or_firm, mycoo, rf_id, 
                         employer_assign,  exec_year, exec_dt)] %>% unique()
```

We didn't use the variable "employer assign" constructed by Marco (2015). Our method gives a more refined sample than does his, and it serves our purpose better.

We generate a random sample from the rejected assignments, ensuring that most filtered assignments are individual-firm assignments.

```{r}
sample_notCOO_transactions = sample_n(assg_or_ee[mycoo==0], 1000) 
View(sample_notCOO_transactions[,.(or_name, ee_name, rf_id)])
```




```{r}
# Construct the set of reassigned patents
assg_doc = left_join(assg_doc, assg_coo)
assg_doc = assg_doc %>% data.table()
assg_doc = assg_doc[order(patent_id, record_dt)]

patent_mycoo = assg_doc[
  order(patent_id, -exec_year)][mycoo==1][ # keep the last transactions of each patent, 
    # it helps to remove employer assignments.
    !duplicated(patent_id),.(patent_id, exec_year, rf_id, mycoo, exec_dt)][
      usfirm_patent[,.(appyear, patent_id, gyear, 
                       enter_year, filedby1perc2, filedby10perc2, 
                       cumN, N5yr, assignee_id)], 
      on = 'patent_id']

patent_mycoo = assignee_cumN[,.(assignee_id, top01allyear, 
                               top1allyear, top10allyear)] %>%
  inner_join(patent_mycoo, .) %>% data.table()

patent_mycoo = assignee_name[patent_mycoo, on = 'assignee_id']
patent_sold = patent_mycoo[mycoo==1]
patent_sold_or_ee = assg_or_ee[,.(ee_dismbg, or_dismbg, rf_id)][patent_sold, on = 'rf_id']

patent_sold_matched = patent_sold_or_ee %>% 
  semi_join(., assignee_name, by = c("ee_dismbg" = 'org_dismbg')) %>%
  data.table()
```

```{r}
nrow(patent_sold_matched[!duplicated(patent_id)])
nrow(patent_sold[!duplicated(patent_id)])
```
575359/794566 patents remained.

```{r}
# calculate the size of buers
cross_table_ee = patent_sold_matched[,.(ee_dismbg, exec_year)] %>% unique() %>% inner_join(usfirm_patent[!is.na(appyear),.(ee_dismbg = org_dismbg, appyear)])
cross_table_ee = data.table(cross_table_ee)
cross_table_ee[,cumN_ee:=sum(exec_year>=appyear), by = .(exec_year, ee_dismbg)]
cross_table_ee[,N5yr_ee:=sum(exec_year>=appyear & appyear>exec_year-5), by = .(exec_year, ee_dismbg)]
cross_table_ee = cross_table_ee %>% select(-appyear) %>% unique()
patent_sold_matched = right_join(cross_table_ee, patent_sold_matched, 
                                 by = c('ee_dismbg', 'exec_year'))
patent_sold_matched = patent_sold_matched%>% data.table()
remove(cross_table_ee)
```

```{r}
ee_ctop = assignee_cumN %>% left_join(., assignee_name) %>% 
  select(ee_dismbg = org_dismbg, ee_top1allyear = top1allyear, 
         ee_organization = organization,
         ee_assignee_id = assignee_id,
         ee_top01allyear = top01allyear, 
         ee_top10allyear = top10allyear) %>% data.table()

ee_ctop[, idN:=.N,by = ee_dismbg]
ee_ctop = ee_ctop[order(ee_dismbg)]
# There are several firms have multiple assignees. By eyebrowsing, these should be the same firms. Let's combine them.
ee_ctop = ee_ctop %>% group_by(ee_dismbg) %>%
  summarize(ee_top01allyear = max(ee_top01allyear),
            ee_top10allyear = max(ee_top10allyear),
            ee_top1allyear = max(ee_top1allyear))

patent_sold_matched = inner_join(patent_sold_matched, ee_ctop)
patent_sold_matched = data.table(patent_sold_matched)

patent_sold_matched = patent_sold_matched %>% 
  inner_join(pat_distr_perfirm[,.(year, q90_5yr, q99_5yr, 
                                  q90_5yr2, q99_5yr2)], by = c('exec_year' = 'year'))
patent_sold_matched = data.table(patent_sold_matched)
patent_sold_matched[, ee_top1:= N5yr_ee > q99_5yr]
patent_sold_matched[, ee_top10:= N5yr_ee > q90_5yr]

patent_sold_matched[, ee_top01allyear:=max(ee_top01allyear), by = rf_id]
patent_sold_matched[, ee_top10allyear:=max(ee_top10allyear), by = rf_id]
patent_sold_matched[, ee_top1allyear:=max(ee_top1allyear), by = rf_id]
patent_sold_matched[, ee_top1:=as.logical(max(ee_top1)), by = rf_id]
patent_sold_matched[, ee_top10:=as.logical(max(ee_top10)), by = rf_id]
patent_sold_matched[,N5yr_ee:=max(N5yr_ee), by = rf_id]

patlevel = patent_sold_matched %>%
  select(contains("ee_top"), contains("N5yr_ee"), rf_id, patent_id, exec_year, appyear, gyear) %>% unique()
patlevel = inner_join(usfirm_patent[,.(gyearN = .N), by = gyear], patlevel) %>%
  inner_join(usfirm_patent[,.(appyearN = .N), by = appyear]) %>%
  inner_join(patent_sold[,.(soldNgyear = .N), by = gyear]) %>%
  data.table()
```


```{r}
patlevel[gyear %in% 1983:2017, mean(ee_top1), by = gyear] %>%
  ggplot(aes(gyear, V1)) + geom_line()+geom_point() + 
  labs(title = "the share of top 1% buyers in matched reassigned patents", x = "Year of application", y = "")
ggsave("top1buer_share.png", width = 10, height = 6)
patlevel[gyear %in% 1983:2017, .(top1buyershare=mean(ee_top1)), by = gyear] %>%
  write.dta(., "dta_for_figures/top1buer_share.dta")
```



```{r}
fee = read_delim("../raw_data/MaintFeeEvents_20180402.txt", delim = ' ',
           col_names = c('patent', 'application', 'small', 'appdate', 'gdate', 'feedate', 'code'),
           col_types = cols_only(patent = 'i', small = 'c', feedate = 'c', code = 'c')) %>% data.table()
fee[, feedate := as.Date(feedate, "%Y%m%d")]
fee = unique(fee)
fee = fee[order(patent, feedate)]
fee[,feeyear:=year(feedate)]
save(fee, file = "../raw_data/fee.RData")
load( "../raw_data/fee.RData")
```

```{r}
big_code = "BIG.|F17[3-5]|F177|LSM[1-3]|M155[1-6]|M18[1-6]|R155[1-6]|R18[1-6]|STOL|STOM"
small_code = "F27[3-5]|F277|M255[1-6]|M27[3-5]|M277|M28[1-6]|M355[1-6]|MICR|R255[1-8]|R27[3-7]|R28[1-6]|R355[1-6]|SM0[1-3]|SMAL"
fee[str_detect(code, big_code), big := 1]
fee[str_detect(code, small_code), big := 0]
```


```{r}
fee_size = fee[!is.na(big)]
fee_size[, patent_id := as.character(patent)]
patent_sold_fee_size = inner_join(patent_sold, fee_size, by = c("patent_id"))
patent_sold_fee_size = data.table(patent_sold_fee_size)
patent_sold_fee_size = patent_sold_fee_size[order(patent_id, exec_dt, feedate)]
patent_sold_fee_size[, Lfeedate :=shift(feedate), by = .(patent_id, exec_dt)]
patent_sold_fee_size[, Lbig:= shift(big), by = .(patent_id, exec_dt)]
```


```{r}
patent_sold_fee_size[, less_feedate:= exec_dt <= feedate]
patent_sold_fee_size[, larger_feedate:= exec_dt > Lfeedate]
patent_sold_fee_size[is.na(larger_feedate), larger_feedate:=T]
patent_sold_fee_size[, or_size_flag := less_feedate & larger_feedate]
patent_sold_or_size = patent_sold_fee_size[or_size_flag==1]

patent_sold_or_size[exec_year >=1992 & exec_year <=2015, sum(1-big)/.N, by = exec_year]%>% ggplot(aes(exec_year, V1)) + geom_line() + geom_point()
patent_sold_or_size[Lbig==0][exec_year >=1992 & exec_year <=2015, sum(1-big)/.N, by = exec_year]%>% ggplot(aes(exec_year, V1)) + geom_line() + geom_point() + labs(x = "Year of reassignment", y = '', title = 'Share of small buyers (conditional on small sellers)')
ggsave("share_of_small_buyers.png", width = 10, height = 6)
patent_sold_or_size[Lbig==0][exec_year >=1992 & exec_year <=2015, sum(1-big)/.N, by = exec_year] %>%
  write.dta(., "dta_for_figures/share_of_small_buyers.dta")
```

```{r}
nrow(patent_sold_or_size[Lbig==0])
```

Only 26287 patents remain.

# Self Citation


```{r}
rm(list = ls())
load(file = "usfirm_patent_1963_2017.RData")
citation = read_delim("../raw_data/uspatentcitation.tsv", delim = "\t",
           col_types = cols_only(patent_id = "c", citation_id = "c"))
save(citation, file = "../raw_data/citation.RData")
# load("../raw_data/citation.RData")
usfirm_patent_1963_2017 = usfirm_patent_1963_2017 %>% select(patent_id, assignee_id, gyear, gdate, appdate, appyear, assignee_type)
setnames(usfirm_patent_1963_2017, old = "patent_id", new = 'patent')
```

```{r}
citation_ij = citation %>% rename(patent_i = patent_id, patent_j = citation_id) %>% 
  mutate(patent_i = str_remove_all(patent_i, '\\s'), patent_j = str_remove_all(patent_j, '\\s'))
  data.table()
# citation_ij = citation_ij[complete.cases(citation_ij)]
colname = colnames(usfirm_patent_1963_2017)
colname_i = sapply(colname, function(x) str_c(x, "_i"))
colnames(usfirm_patent_1963_2017) = colname_i
citation_ij = inner_join(usfirm_patent_1963_2017, citation_ij)
# Use inner_join, only a 1% was deleted from citation.

colname_j = sapply(colname, function(x) str_c(x, "_j"))
colnames(usfirm_patent_1963_2017) = colname_j
citation_ij = left_join(citation_ij, usfirm_patent_1963_2017)
citation_ij = citation_ij %>% data.table()
colnames(usfirm_patent_1963_2017) = colname
```

```{r}
citation_usfi = citation_ij[assignee_type_i == '2']
citation_usfi = citation_usfi[appyear_i >=1975 & appyear_i<=2014]
citation_usfij = citation_usfi[assignee_type_j == '2']
citation_usfij5yr = citation_usfi[appyear_i <= appyear_j+5]
citation_usfij10yr = citation_usfi[appyear_i <= appyear_j+10]
```

Be careful of potential selection issue:

```{r}
usfi_pat_n = citation_usfi[!duplicated(patent_i), .N, by = appyear_i]
usfij5yr_pat_n = citation_usfij5yr[!duplicated(patent_i), .(N5yr = .N), by = appyear_i]
usfij10yr_pat_n = citation_usfij10yr[!duplicated(patent_i), .(N10yr = .N), by = appyear_i]
inner_join(usfi_pat_n, usfij5yr_pat_n) %>% inner_join(usfij10yr_pat_n) %>%
  ggplot(aes(appyear_i)) + geom_line(aes(y = N, color = 'Full')) +
  geom_line(aes(y = N5yr, color = '5yr')) + 
  geom_line(aes(y = N10yr, color = '10yr'))
inner_join(usfi_pat_n, usfij5yr_pat_n) %>% inner_join(usfij10yr_pat_n) %>% 
  write.dta("dta_for_figures/sample_size_citation_data.dta")
```


```{r}
citation_usfi[,self_cite := assignee_id_i == assignee_id_j]
citation_usfi[appyear_i>=1980, sum(self_cite, na.rm = T)/.N, by = appyear_i] %>%
  ggplot(aes(appyear_i, V1))+ geom_line() + geom_point() + labs(title = "Full Sample: Self Citation Share")
ggsave("full_sample_self_citation.png",width = 10, height = 6)
citation_usfi[appyear_i>=1980, .(self_cite_share = sum(self_cite, na.rm = T)/.N), by = appyear_i] %>%
  write.dta(.,'dta_for_figures/full_sample_self_citation.dta')

citation_usfij5yr[,self_cite := assignee_id_i == assignee_id_j]
citation_usfij5yr[appyear_i>=1980, sum(self_cite, na.rm = T)/.N, by = appyear_i] %>%
  ggplot(aes(appyear_i, V1))+ geom_line() + geom_point() + labs(title = "5yr Window: Self Citation Share")
ggsave("5yr_window_self_citation.png",width = 10, height = 6)
citation_usfij5yr[appyear_i>=1980, .(self_cite_share = sum(self_cite, na.rm = T)/.N), by = appyear_i] %>%
  write.dta(.,'dta_for_figures/5yr_window_self_citation.dta')

citation_usfij10yr[,self_cite := assignee_id_i == assignee_id_j]
citation_usfij10yr[appyear_i>=1980, sum(self_cite, na.rm = T)/.N, by = appyear_i] %>%
  ggplot(aes(appyear_i, V1))+ geom_line() + geom_point() + labs(title = "10yr Window: Self Citation Share")
ggsave("10yr_window_self_citation.png",width = 10, height = 6)
citation_usfij10yr[appyear_i>=1980, .(self_cite_share = sum(self_cite, na.rm = T)/.N), by = appyear_i] %>%
  write.dta(.,'dta_for_figures/10yr_window_self_citation.dta')
```



# Claim Length

```{r}
patent_claim = read.dta("../raw_data/patent-claim-length.dta") %>% data.table()
patent_claim[, patent_id :=  str_remove_all(patnum, '\\s')]
load("usfirm_patent.RData")
patent_claim = patent_claim %>% inner_join(usfirm_patent) %>% data.table()
```


```{r}
patent_claim[appyear>=1975 & appyear <=2014, .(mean(claim_length, na.rm = T)), by = appyear] %>% ggplot(aes(appyear)) +
  geom_line(aes(y = V1)) + geom_point(aes(y = V1)) + labs(y = '', x = 'Year of application',  title = "Average Claim Length")
ggsave("mean_claim_length.png", width = 10, height = 6)

patent_claim[appyear>=1975 & appyear <=2014, .(mean_claim_length = mean(claim_length, na.rm = T)), by = appyear] %>% 
  write.dta(., "dta_for_figures/mean_claim_length.dta")

```


