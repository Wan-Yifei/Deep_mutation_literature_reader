import csv
import math
import sys
from Bio import Entrez
from pickle_toolbox import *


def get_pmid(input_file):
    """
    Extract pmids from input table file.

    :param input_file: csv file has pmids in header;
    :return: a list of pmid.
    """
    pmids = []
    with open(input_file, "r", encoding='utf-8') as csvfile:
        raw = csv.reader(csvfile, delimiter=",")
        line = next(raw)
        pmid_index = line.index("pmids")
        line = next(raw)
        while line:
            pmid = line[pmid_index]
            # print(pmid)
            if pmid:
                if "|" in pmid:
                    pmids += pmid.split("|")
                else:
                    pmids.append(pmid)
            try:
                line = next(raw)
            except StopIteration:
                line = ""
        pmids = list(set(pmids))  # remove duplicates
    return pmids


def get_pickled_pmid(input_file, hmac_key, use_hmac):
    if use_hmac:
        pmids_raw = pickle_load(input_file, hmac_key=hmac_key, use_hmac=use_hmac)
    else:
        pmids_raw = pickle_load(input_file)
    if isinstance(pmids_raw, list):
        pmids = pmids_raw
        return pmids
    if isinstance(pmids_raw, dict):
        pmids = [pmid for pmid in pmids_raw.keys()]
        return pmids


def fetch_title_abstract(pmids, email, tool_nm="fetch_abstract_api", api_key="", max_tries=3, sleep_between_tries=15,
                         pmids_limit=200):
    """
    Fetch title and abstract of each PMID. If PMID in the response does not match the PMID in the query would be record
    as mismatch pmid.

    :param pmids: The list of all PMIDs;
    :param email: The e-mail addresss of user for NCBI API access;
    :param tool_nm: Tool name of the tool (current program) for NCBI API access;
    :param api_key: Personal API key from NCBI. If not set, 3 queries/second are allowed. 10 queries/second with a valid
            API key;
    :param max_tries: How many times failed requests will be automatically retried on error;
    :param sleep_between_tries: The delay, in seconds, before retrying a request on error;
    :param pmids_limit: the limit of the number of PMIDs per each query;
    :return: pmids_abstract: dictionary key -> pmid, value -> abstract;
    :return: pmids_title: dictionary key -> pmid, value -> title;
    :return: mismatch_pmid: list of mismatch PMIDs.
    """
    pmids_split = [pmids[ind * pmids_limit: (ind + 1) * pmids_limit] for ind in
                   range(math.ceil(len(pmids) / pmids_limit))]  # Every query has pmids_limits # pmids.
    pmids_abstract = {}
    pmids_title = {}
    mismatch_pmids = []
    # Config Entrez
    Entrez.email = email
    Entrez.tool = tool_nm
    Entrez.api_key = api_key
    Entrez.max_tries = max_tries
    Entrez.sleep_between_tries = sleep_between_tries
    # Fetch abstract
    batch = 0  # batch count of query, per pmids_limts pmis / one batch.
    for pmid_list in pmids_split:
        batch += 1
        if len(pmid_list) > 1:
            pmid_list = ",".join(pmid_list)  # Use commas to split multiple pmids
        else:
            pmid_list = pmid_list[0]  # flat the list to get single pmid
        handle = Entrez.efetch(db="pubmed", retmode="xml", id=pmid_list)
        records = Entrez.read(handle)
        response_pmids = [str(records["PubmedArticle"][ind]["MedlineCitation"]["PMID"]) for ind
                          in range(len(records["PubmedArticle"]))]
        if set(response_pmids) == set(pmid_list.split(",")):
            print("Response {} PMIDs in the batch {} match submitted PMIDs, QC pass!\n".format(len(set(response_pmids)),
                                                                                               batch))
        else:
            # record mismatch PMIDs for manual check!
            mismatch_pmids += [pmid for pmid in set(pmid_list.split(",")) if pmid not in response_pmids]
            print("Response {} PMIDs in the batch {} do not match submitted PMIDs, QC failed\n"
                  .format(len(set(response_pmids)), batch))

        for record in records["PubmedArticle"]:
            if str(record["MedlineCitation"]["PMID"]) not in mismatch_pmids:
                try:
                    pmids_abstract[str(record["MedlineCitation"]["PMID"])] = str(record["MedlineCitation"]["Article"]
                                                                                 ["Abstract"]["AbstractText"][0])

                    pmids_title[str(record["MedlineCitation"]["PMID"])] = str(record["MedlineCitation"]["Article"]
                                                                              ["ArticleTitle"])
                except KeyError:
                    pmids_abstract[str(record["MedlineCitation"]["PMID"])] = ""
                    pmids_title[str(record["MedlineCitation"]["PMID"])] = ""
        handle.close()
    for pmid in mismatch_pmids:
        pmids_abstract[pmid] = ""
        pmids_title[pmid] = ""
    print("Please check mismatch PMIDs: {}".format(mismatch_pmids))
    return pmids_abstract, pmids_title, mismatch_pmids


def abstract_status_check(pmids_abstract, mismatch_pmids):
    """
    Check status of abstract of each PMID.

    :param pmids_abstract: dictionary key -> pmid, value -> abstract;
    :param mismatch_pmids: list of mismatch PMIDs;
    :return: abstract_status: dictionary key -> pmid, value -> status.
    """
    abstract_status = {}
    for pmid in pmids_abstract.keys():
        if pmid in mismatch_pmids:
            abstract_status[pmid] = "PMID mismatch"
        elif pmids_abstract[pmid] == "":
            abstract_status[pmid] = "Abstract not found"  # Not valid abstract in the response
        else:
            abstract_status[pmid] = "Abstract ready"
    return abstract_status


def title_status_check(pmids_title, mismatch_pmids):
    """
    Check status of title of each PMID.

    :param pmids_title: dictionary key -> pmid, value -> title;
    :param mismatch_pmids: list of mismatch PMIDs;
    :return: title_status: dictionary key -> pmid, value -> status.
    """
    title_status = {}
    for pmid in pmids_title.keys():
        if pmid in mismatch_pmids:
            title_status[pmid] = "PMID mismatch"
        elif pmids_title[pmid] == "":
            title_status[pmid] = "Title not found"  # Not valid title in the response
        else:
            title_status[pmid] = "Title ready"
    return title_status


def summarize_status(dict_status, status_type):
    """
    Summrize and print out status check result of title or abstract.

    :param dict_status: dictionary key -> pmid, value -> title or abstract;
    :param status_type: specify the type of input dictionary, "abstract" or "title";
    :return: None
    """
    count_ok = 0
    count_mismatch = 0
    count_no_found = 0
    if "abstract" == status_type:
        status_type = "abstract"
        for pmid in dict_status.keys():
            if dict_status[pmid] == "PMID mismatch":
                count_mismatch += 1
            if dict_status[pmid] == "Abstract not found":
                count_no_found += 1
        count_ok = len(dict_status) - count_mismatch - count_no_found
    elif "title" == status_type:
        status_type = "title"
        for pmid in dict_status.keys():
            if dict_status[pmid] == "PMID mismatch":
                count_mismatch += 1
            if dict_status[pmid] == "Abstract not found":
                count_no_found += 1
        count_ok = len(dict_status) - count_mismatch - count_no_found
    else:
        status_type = None
        print("Error: Status type must be named as abstract or title!")
    if status_type:
        print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("Summary of {} status:".format(status_type))
        print("Fetch {} successfully: {}".format(status_type, count_ok))
        print("PMID mismatch count: {}".format(count_mismatch))
        print("Not found {}: {}".format(status_type, count_no_found))
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")


def pickle_output(content, output_file, output_folder, key="", use_hmac=False):
    """
    Pickle content as bytes file to specified location.

    :param content: object would be pickled;
    :param output_file: output file name;
    :param output_folder: output folder, would be the root of output path;
    :param key: key for hash message authentication;
    :param use_hmac: Apply hash message authentication to pickled file or not, True or False. False by default;
    :return: None
    """
    output_path = output_folder + "/" + output_file
    pickled_content = pickle.dumps(content)
    if use_hmac:
        hmac_header = generate_hmac_signature(pickled_content, key)
        print("Hash message signature: {}".format(hmac_header))
        hmac_header = hmac_header + "\n"  # concatenate HMAC signature to a new line
    else:
        hmac_header = None
    with open(output_path, "w+") as output:
        if hmac_header:
            output.write(hmac_header)
    with open(output_path, "ab") as output:
        output.write(pickled_content)
        print("\nPickled file to {}".format(output_path))


def main():
    input_file = sys.argv[1]
    output_folder = sys.argv[2]
    email = sys.argv[3]
    try:
        hmac_key = sys.argv[4]
    except IndexError:
        hmac_key = None
    if hmac_key:
        use_hmac = True  # When HMAC key is provided, use HMAC to sign pickled output.
    else:
        use_hmac = False
    if input_file.endswith("csv"):
        pmids = get_pmid(input_file)
    if input_file.endswith("pickle"):
        pmids = get_pickled_pmid(input_file, hmac_key, use_hmac)
    pmids_abstract, pmids_title, mismatch_pmids = fetch_title_abstract(pmids, email)
    abstract_status = abstract_status_check(pmids_abstract, mismatch_pmids)
    title_status = title_status_check(pmids_title, mismatch_pmids)
    summarize_status(abstract_status, status_type="abstract")
    summarize_status(title_status, status_type="title")
    print("\nSave PMID-TITLE dictionary:")
    pickle_output(pmids_title, "pmid_title.pickle", output_folder, key=hmac_key, use_hmac=use_hmac)
    print("\nSave PMID-ABSTRACT dictionary:")
    pickle_output(pmids_abstract, "pmid_abstract.pickle", output_folder, key=hmac_key, use_hmac=use_hmac)


if __name__ == "__main__":
    main()
