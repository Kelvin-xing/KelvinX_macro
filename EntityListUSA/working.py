### This code cleans Code of Federal Regulations Supplement No. 4 to Part 744, Title 15, from https://www.ecfr.gov/current/title-15/subtitle-B/chapter-VII/subchapter-C/part-744/appendix-Supplement%20No.%204%20to%20Part%20744
### It contains the list of names of certain non-US entities that are subject to specific license requirements for export, re-export or transfer of specified items by bureau of industry and security, US Department of Commerce
### CSV version downloadable from https://www.bis.doc.gov/index.php/documents/consolidated-entity-list?format=html
### For more details, check https://www.bis.doc.gov/index.php/policy-guidance/lists-of-parties-of-concern/entity-list
### A search engine is available on https://www.trade.gov/consolidated-screening-list

import pandas as pd
import re
from dateutil.parser import parse
from collections import OrderedDict


## Historical revisions updated 2023-12-22
dates = ['2023-12-07','2023-11-21','2023-11-17','2023-11-06','2023-10-19','2023-10-11','2023-9-27','2023-7-19','2023-6-21','2023-6-14','2023-5-22','2023-4-26','2023-4-17','2023-3-30','2023-3-14','2023-3-06','2023-2-27','2023-2-14','2023-2-01','2022-12-23','2022-12-19','2022-12-16','2022-12-08','2022-10-21','2022-10-13','2022-10-07','2022-10-04','2022-9-16','2022-9-09','2022-8-24','2022-6-30','2022-6-06','2022-6-01','2022-5-11','2022-4-11','2022-4-07','2022-3-16','2022-3-09','2022-3-08','2022-3-03','2022-2-14','2022-2-03','2021-12-17','2021-11-26','2021-11-04','2021-10-05','2021-8-20','2021-7-21','2021-7-19','2021-7-12','2021-7-06','2021-6-24','2021-6-16','2021-6-01','2021-4-09','2021-3-16','2021-3-08','2021-3-04','2021-1-15','2020-12-23','2020-12-22','2020-10-30','2020-10-19','2020-9-22','2020-9-11','2020-8-27','2020-8-20','2020-7-22','2020-6-18','2020-6-05','2020-5-19','2020-3-16','2020-3-02','2019-12-18','2019-12-06','2019-11-13','2019-10-21','2019-10-09','2019-8-21','2019-8-14','2019-6-24','2019-5-24','2019-5-21','2019-5-14','2019-4-11','2018-12-20','2018-10-30','2018-9-26','2018-9-12','2018-9-04','2018-8-30','2018-8-01','2018-3-22','2018-2-16','2018-1-26','2017-12-20','2017-9-25','2017-6-30','2017-6-22','2017-5-26','2017-4-18','2017-3-29','2017-3-16']

## File under cleaning process
# Note: Can loop over dates
file = f"EntityList{dates[0]}.csv"
df = pd.read_csv(f"PolicyData/EntityListUSA/{file}",usecols=['Name', 'Address','City','State/Province','Country','Federal Register Notice','Effective Date','License Requirement','License Policy','Alternate Name']).astype(str).replace('nan','')

## Join identifier
# Note: the identifier is from Consolidated Screening List from US government's International Trade Administration on https://www.trade.gov/consolidated-screening-list
# Note: the identifier is generally unique to name, but there are still entities with multiple addresses assigned different identifiers. I have manually checked and it appears to be the dataset's problem and they essentially refer to the same entities. Duplicate identifiers are dropped and only the first are kept.
idmatch = (pd.read_csv(f"PolicyData/EntityListUSA/consolidated2024-01-05.csv",usecols=['_id','name'])
            .astype(str)
            .drop_duplicates(subset="name",keep='first')
            .set_index('name'))
df = df.join(idmatch, on='Name', how='inner')

## Cleaning Name
# standardize
replacer = OrderedDict({
            '?':'',
            "Xi'an":'Xian',
            ", a.k.a, the following ten aliases":'',
            ', a.k.a.':'; ',
            ', and subordinate entity Nuclear reactors (including power plants), fuel reprocessing and enrichment facilities, all uranium processing, conversion and enrichment facilities, heavy water production facilities and any collocated ammonia plants.':'',
            ', and subordinate entity':';',
            ' (Parent Organization: China National Nuclear Group Corporation (CNNC))':'',
            '¬†':' ',
            ', including the Armed Forces of Belarus and all operating units wherever located.  This includes the national armed services (army and air force), as well as the national guard and national police, government intelligence or reconnaissance organizations of the Republic of Belarus.  All addresses located in Belarus.':'',
            ' (Parent Organization: China National Nuclear Group Corporation (CNNC))':'',
            ' ‚ÄúState Machine Building Design Bureau ‚ÄúVympel‚Äù By Name I.I.Toropov‚Äù':'',
            ' (f.k.a., Federalnoe Gosudarstvennoe Byudzhetnoe Uchrezhdenie Sanatori Nizhnyaya Oreanda Upravleniya)':'',
            '‚Äú':' ',
            '‚Äù':' ',
            '‚Äì':' ',
            '            and , Sharjah, U.A.E.':''})
# to conduct the standardization
def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text
df['Name'] = df['Name'].str.strip().apply(lambda x:replace_all(x,replacer))
# concatenate with alternate names to expand into panel later
df['Names'] = df.apply(lambda row: row['Name'] + '; ' + row['Alternate Name'] if row['Alternate Name'] else row['Name'], axis=1)
# drop alternate name
df = df.drop(['Alternate Name','Name'],axis=1)


## Concatenate addresses
df['Full Address'] = (df['Address'] + ' '+ df['City'] + ' ' + df['State/Province']).str.strip()
df = df.drop(['Address','City','State/Province'],axis=1)


## Explode into panel
df = df.set_index(['Names','Full Address','Country','License Requirement','License Policy'])
# df = df.set_index(['Effective Date','Federal Register Notice','Full Address','Country','License Requirement','License Policy'])
df['count_date'] = df['count_register'] = 0
# splitting
for index in range(len(df)):
    df['Effective Date'][index] = re.split(r'[;,]\s*',df['Effective Date'][index])
    df['count_date'][index] = len(df['Effective Date'][index])
    df['Federal Register Notice'][index] = re.split(r'\sand\s|[;,]\s*',df['Federal Register Notice'][index])
    df['count_register'][index] = len(df['Federal Register Notice'][index])
    # df['Names'][index] = re.split(r';and\s*|;\s*',df['Names'][index])
# exploding
df = (df
        # .explode('Effective Date')
        # .explode('Federal Register Notice')
        # .explode('Names')
        .reset_index())


## Drop duplicates and re-order
df = (df
    #   .reindex(columns=['_id','Names','Full Address','Type','Country','Federal Register Notice','Effective Date','License Requirement','License Policy'])
      .reindex(columns=['_id','Names','Full Address','Type','Country','count_register','count_date','License Requirement','License Policy'])
      .drop_duplicates()
      .reset_index(drop=True))

## Save to a new csv file
file_cleaned = file.replace('.csv','_cleaned_test.csv')
df.to_csv(f"PolicyData/EntityListUSA/{file_cleaned}")