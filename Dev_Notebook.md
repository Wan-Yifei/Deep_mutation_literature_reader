# Dev Notebook

## Resources:

1. Biomedical natural language processing\
http://bio.nlplab.org/
   
2. BioWordVec & BioSentVec: pre-trained embeddings for biomedical words and sentences\
https://github.com/ncbi-nlp/BioSentVec
   
3. PDFminer: \
   https://pdfminer-docs.readthedocs.io/programming.html
   
4. New MVL data: \
data/raw_mvl_data.csv
   

## To do list:

1. ~~need to convert PMID to DOI~~ -> Done

Use Bio.Entrez to search dois for pmids. But a part of pmids don't have
corresponding dois. Set doi as "" when a pmid does not have
valid doi.

2. ~~download paper from Sci-Hub~~ -> Done


reference repo:
https://github.com/gadilashashank/Sci-Hub/blob/master/sci_hub.py

+ download summary:

```
20121/08/30
Summary of PMID-DOI source file:

Total # pmids: 3598\
Total # doi: 3467

Count # of ok: 3330\
Count # of captcha encountered: 143\
Count # of status unknown: 125
```
2.1 Only download abstract
```
2021/11/17 
source: raw_mvl_data.csv
pmid file: pmid_source_new.pickle
mismatch PMIDs: ['30297953', '33046410', '187169170', '92582499', '2249', '20301355', '245490550', '274435140', '20301779', '20301425', '22720333', '20301428', '121769658', '23576526', '23035301', '22379635', '33119245', '27963154', '27963596', '20301519', '22377370', '20301572', '20301348', '20301679
', '20301575', '28022460', '184150', '1092248265']

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Summary of abstract status:
Fetch abstract successfully: 15392
PMID mismatch count: 28
Not found abstract: 612
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Summary of title status:
Fetch title successfully: 16004
PMID mismatch count: 28
Not found title: 0
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

```

3. ~~Limit the request times per second~~ -> Done

Both NCBI query and Sci-Hub download should be limited.
+ ~~Update Bio to v1.7.8 to automatically limit the frequency of DOI queries.~~
+ ~~Manually limit the download by `time.sleep`.~~

4. When Sci-Hub does not provide required paper, may try to fetch the
paper from Pubmed directly?


5. Convert PDF to text\
+ convert full text;
+ extract title and abstract;
+ extract body.

6. ~~Get abstract text and title from Entrez~~ -> Done
+ ~~get abstract~~
+ ~~Get title~~ 
+ ~~pickle output~~
+ ~~write output to txt~~
+ ~~Apply HMAC to pickle~~

7. Data clean
+ ~~remove copyright symbol;~~
+ ~~capture HGVS variants by regex before removing punctuation.~~
+ integrate NER tool tmVar to extract variant from literature\
https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmids=24451227&concepts=mutation
  
## Test cmd

+ Extract PMID-DOI

E:\anaconda\python.exe extract_pmid_doi.py data\MVL_HC_ongoing.csv data\pmid_doi_list.csv wanyifei0123@gmail.com

+ Download PDF

E:\anaconda\python.exe download_pmid.py data\pmid_doi_list.csv data\pmid_refs

+ Get abstract

E:\anaconda\envs\text_mining\python.exe E:/project/Deep_mutation_literature_reader/fetch_pmid_title_abstract.py E:\project\Deep_mutation_literature_reader\data\MVL_HC_ongoing.csv .\data\title_abstract_refs wanyifei0123@gmail.com wanyifei


## Issue list:

1. some PMIDs cannot be found in PUBMED.
2. ~~Doi not match title content. Title from `soup` represents `>` as `&gt` somewhere.~~

Fixed by `soup.find("").encode(formatter=None)`.

e.g.

URL `https://sci-hub.se/10.1002/1098-1004(200011)16:5%3C417::AID-HUMU6%3E3.0.CO;2-4`

DOI `10.1002/1098-1004(200011)16:5<417::AID-HUMU6>3.0.CO;2-4`

Title `<title>Sci-Hub | Mutation analysis of the hamartin gene using denaturing high performance liquid chromatography. Human Mutation, 16(5), 417–421 | 10.1002/1098-1004(200011)16:5&lt;417::AID-HUMU6&gt;3.0.CO;2-4</title>
10.1002/1098-1004(200011)16:5<417::AID-HUMU6>3.0.CO;2-4`

3. ~~Sci-Hub captcha encountered~~

Both get url and download pdf have this issue.\
Fix: skip the pmid while captcha encountered.

4. Data conflict:
One pmid is pathogenic and benign at the same time.
   e.g. `24728327`
   In MVL_HC_ongoing.csv line: 239