import sys
import pickle
import contractions
from fetch_pmid_title_abstract import pickle_output
from hgvs_variants_capture import *
from data_clean_toolkit import *
from hmac_signature import *


def main():
    input_title = sys.argv[1]
    input_abstract = sys.argv[2]
    hmac_key = sys.argv[3]
    output_folder = sys.argv[4]
    with open(input_title, "rb") as title:
        title_signature_old = title.readline().decode("utf-8")
        title_content = b"".join(title.readlines())
        if verify_hmac_signature(title_content, hmac_key, title_signature_old):
            pmids_title = pickle.loads(title_content)
            print("Dict pmids_title is ready!")
    with open(input_abstract, "rb") as abstract:
        abstract_signature_old = abstract.readline().decode("utf-8")
        abstract_content = b"".join(abstract.readlines())
        if verify_hmac_signature(abstract_content, hmac_key, abstract_signature_old):
            pmids_abstract = pickle.loads(abstract_content)
            print("Dict pmids_abstract is ready!")
    # 1. Lower all cases
    pmids_title_l = {pmid: title.lower() for pmid, title in pmids_title.items()}
    pmids_abstract_l = {pmid: abstract.lower() for pmid, abstract in pmids_abstract.items()}
    # 2. Expand the contractions
    pmids_title_lc = {pmid: contractions.fix(title) for pmid, title in pmids_title_l.items()}
    pmids_abstract_lc = {pmid: contractions.fix(abstract) for pmid, abstract in pmids_abstract_l.items()}
    # 3. Remove URL
    pmids_title_lcu = {pmid: remove_url(title) for pmid, title in pmids_title_lc.items()}
    pmids_abstract_lcu = {pmid: remove_url(abstract) for pmid, abstract in pmids_abstract_lc.items()}
    # 4. Remove HTML
    pmids_title_lcuh = {pmid: remove_html(title) for pmid, title in pmids_title_lcu.items()}
    pmids_abstract_lcuh = {pmid: remove_html(abstract) for pmid, abstract in pmids_abstract_lcu.items()}
    # 5. Remove special characters
    pmids_title_lcuhs = {pmid: remove_special_characters(title) for pmid, title in pmids_title_lcuh.items()}
    pmids_abstract_lcuhs = {pmid: remove_special_characters(abstract) for pmid, abstract in pmids_abstract_lcuh.items()}
    # 6. Sentence level tokenize
    pmids_title_lcuhst = {pmid: title.split(". ") for pmid, title in pmids_title_lcuhs.items()}
    pmids_abstract_lcuhst = {pmid: abstract.split(". ") for pmid, abstract in pmids_abstract_lcuhs.items()}
    # 7. Remove punctuations
    vartype = "dna"
    pmids_title_lcuhstp = {pmid: [variant_safe_punctuation_remover(pmid, sen, vartype) for sen in title_sens] for
                           pmid, title_sens in pmids_title_lcuhst.items()}
    pmids_abstract_lcuhstp = {pmid: [variant_safe_punctuation_remover(pmid, sen, vartype) for sen in abstract_sens] for
                              pmid, abstract_sens in pmids_abstract_lcuhst.items()}
    # 8. Output cleaned and tokenized text
    pickle_output(pmids_title_lcuhstp, "cleaned_title.pickle", output_folder, key=hmac_key, use_hmac=True)
    pickle_output(pmids_abstract_lcuhstp, "cleaned_abstract.pickle", output_folder, key=hmac_key, use_hmac=True)


if __name__ == "__main__":
    main()
