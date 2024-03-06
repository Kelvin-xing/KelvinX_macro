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
            '            and , Sharjah, U.A.E.':'',
            '√™':'E',
            '‚Äö':'e',
            '‚Äì':''})
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

## Fill country blank
for index in range(len(df)):
    if not df['Country'][index]:
        df['Country'][index] = 1# retrieve the final part of full address
    if df['Country'][index] == 'Kowloon':
        df['Country'][index] = 'Hong Kong'
    elif df['Country'][index] == 330002:
        df['Country'][index] = 'China'


## Explode into panel
# df = df.set_index(['Names','Full Address','Country','License Requirement','License Policy'])
df = df.set_index(['Effective Date','Federal Register Notice','Full Address','Country','License Requirement','License Policy'])
# splitting
for index in range(len(df)):
    # df['Effective Date'][index] = re.split(r'[;,]\s*',df['Effective Date'][index])
    # df['Federal Register Notice'][index] = re.split(r'\sand\s|[;,]\s*',df['Federal Register Notice'][index])
    df['Names'][index] = re.split(r';and\s*|;\s*',df['Names'][index])
# exploding
df = (df
        # .explode('Effective Date')
        # .explode('Federal Register Notice')
        .explode('Names')
        .reset_index())



## Split human and company with identifier as a new column
# keywords for firm, government and research institutes
firmwords = [# Typical firm suffix
            'llc','corp','system','industr','company','ltd','co.','group','factory','enterprise','association','jsc','plant','branch','limited','llp','associates','foundation','inc','sdn bhd','development','headquarter','gmbh','limited','private','oao','ooo','zao','s.a.',' ao',' ab',' oy',' sarl','s.a.l.','trust','fze','fzco','holding','sdn','contracting','complex',' ag','Aktsionernoe Obshchestvo','venture',
            # Geography
            'global','beijing','international',
            # Industry
                # construction and manufacture
                'construction', 'steel','engineer','metro','bridge','production','konstrukt','manufactur','fku uprdor','establishment','Al-Qertas','Vangurd Tec','ELPROM','Elara',
                # energy
                'energy','dietsmannnile','nyakek and sons','Oranto Petroleum','sanco',
                # pharma
                'pharm','ELEMED','medical','Xinjiang Silk Road BGI','rau',
                # IT
                'comput','semiconductor','electron','micro','radio','microwave','cloud','display','elec','infotec','video','tronic','chip','cyber','angstrem','interscan','Proven Glory','higon','sugon','IFLYTEK','netposa','sensenets','network','dji','kindroid','candiru','corad','MCST Lebedev','NPP Istok','Avanlane','Milur SA','elektronika','SMT-iLogic','streloy','Grant Instrument','elektro',
                # tech in general
                'technolog','integra','tech','solution','tekno','armyfly',
                # trade
                'trad','service','export','import','logist',
                # transportation
                'aero','airline','aerospace','shipyard','ship','aircraft','flight','aviation','used car','motors','skylink','vehicle','bike','concord','UEC-Saturn','MPI VOLNA','FASTAIR','Aviazapchast PLC','PT Air',
                # chemistry
                'chemie','interlab','labinvest','femteco',
                # others
                'design','field','resort','focus middle east','cosmos','NM-Tekh','TROJANS','Pearl Coral 1173 CC','consulting',
             # Firm specific
             'huawei', 'Moselectronproekt','nexus','proven honour','gazprom','oceanos','vad, ao','proexcom','magnetar','apex','melkom','abris','dm link','sngb ao','jadeshine','sputnik','ikco','cytrox','serop','aviton','bitreit','mekom','mces','meo','satco','chernomorneftegaz','aquanika','stroygazmontazh','transoil','intelcom','zener','source','rosneft','surgutneftegas','yaltinskaya kinodstudiya','otkrytoe aktsionernoe obshchestvo vneshneekonomicheskoe obedinenie tekhnopromeksport','zte','micado','ar kompozit kimya','Regionsnab','Adimir OU','UAB Pella-Fjord','Alfakomponent','The Mother Ark.','Hasa Nederland B.V.','AST Components',
             # AO prefix 
             'AO Kronshtadt','AO Rubin','AO Aviaagregat','AO PKK Milandr','AO Papilon','AO Geomir','AO SET-1',
             ]
institutewords = ['university','research center','institute','academy','laborator','nscc','advanced research','development center','TsKB MT Rubin','SDB IRE RAS']
govwords =['ministry','desto','paec','cnsim','state ','federal','bureau','intelligence','department','committee','defense','atomic',"people's republic",'glavgosekspertiza rossii','crimea','zorsecurity','police']
words_list = [firmwords,institutewords,govwords]
(list(map(str.lower, words)) for words in words_list)

# matching type
# Note: in this order to capture gov or institutes wrongly taken as firm
df['Type']='Human'
for index, row in df.iterrows():
    text = row['Names']
    if any (firmword in text.lower() for firmword in firmwords):
        df.at[index,'Type']='Firm'
    if any (govword in text.lower() for govword in govwords):
        df.at[index,'Type']='Government'
    if any (instituteword in text.lower() for instituteword in institutewords):
        df.at[index,'Type']='Research Institute'
# idiosyncratic tuning
firm_idio = ['The Jordanian Lebanese Company for Laboratory Instruments S.A.L.','Svyaz Design Bureau, OJSC']
df.loc[df['Names'].isin(firm_idio),'Type'] = 'Firm'


## Clean fed registration
##### to be drawn from the main table and match by name
fc = pd.read_excel('PolicyData/EntityListUSA/UnknownFedRegister.xlsx') # fc for fed register check
# function to tell date aside from register
def split_and_parse(string):
    halves = string.split('. ')
    if len(halves) == 2:
        return halves[1]
    else:
        return halves[0]

# df['Federal Register Notice'] = df['Federal Register Notice'].apply(lambda x: split_and_parse(x))

# for i in range(len(df)):
#     # to save only the register part
#     df['Federal Register Notice'][i] = df['Federal Register Notice'][i].split(' No.',1
#     )[0].split(' (',1)[0].replace(' no. 242 pg.','').strip()
#     # maunally fill na
#     if df['Federal Register Notice'][i] == 'nan':
#         continue
#         # match from fc
#     else:
#         continue


## Match Effective Date
# Fed Register - Effective date 2011-2023
# Note: companies before 2011 left unchanged
dc = pd.read_excel('PolicyData/EntityListUSA/FedRegisterDates2011-2023.xlsx','Table_cleaned') # dc short for date check

# standardize date




## Clean date
# date = parse(str(df['Effective Date']))
# df['Effective Date'] = date.strftime("%Y-%m-%d")


## Drop duplicates and re-order
df = (df
      .reindex(columns=['_id','Names','Full Address','Type','Country','Federal Register Notice','Effective Date','License Requirement','License Policy'])
      .drop_duplicates()
      .reset_index(drop=True))

## Save to a new csv file
file_cleaned = file.replace('.csv','_cleaned.csv')
df.to_csv(f"PolicyData/EntityListUSA/{file_cleaned}")