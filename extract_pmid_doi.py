import sys
import csv
import math
from Bio import Entrez  # Bio.Entrez at leagst v1.75 to enable auto rate limit


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
            except:
                line = ""
        pmids = list(set(pmids))  # remove duplicates
    return pmids


def get_doi_api(pmids, email, tool_nm="get_doi_api", api_key="", max_tries=3, sleep_between_tries=15, pmids_limit=200):
    """
    Query DOI for each PMID.
    :param pmids_limit:
    :param pmids: The list of all PMIDs;
    :param email: The e-mail addresss of user for NCBI API access;
    :param tool_nm: Tool name of the tool (current program) for NCBI API access;
    :param api_key: Personal API key from NCBI. If not set, 3 queries/second are allowed. 10 queries/second with a valid
            API key;
    :param max_tries: How many times failed requests will be automatically retried on error;
    :param sleep_between_tries: The delay, in seconds, before retrying a request on error;
    :return: pmids_doi: a dictionary with key -> PMID and value -> DOI..
    """
    pmids_split = [pmids[ind * pmids_limit: (ind + 1) * pmids_limit] for ind in
                   range(math.ceil(len(pmids) / pmids_limit))]
    pmids_doi = {}
    # Config Entrez
    Entrez.email = email
    Entrez.tool = tool_nm
    Entrez.api_key = api_key
    Entrez.max_tries = max_tries
    Entrez.sleep_between_tries = sleep_between_tries
    invalid_pmid = 0  # count PMIDs cannot be found in PUBMED
    for pmid_list in pmids_split:
        if len(pmid_list) > 1:
            pmid_list = ",".join(pmid_list)  # Use commas to split multiple pmids
        else:
            pmid_list = pmid_list[0]  # flat the list to get single pmid
        handle = Entrez.esummary(db="pubmed", retmode="xml", id=pmid_list)
        try:
            records = Entrez.read(handle)
            for record in records:
                if record["Id"] == record["ArticleIds"]["pubmed"][0]:  # Check the response ID matches input pmid
                    try:
                        pmids_doi[record["ArticleIds"]["pubmed"][0]] = record["ArticleIds"]["doi"]
                    except:
                        pmids_doi[record["ArticleIds"]["pubmed"][0]] = ""
                else:
                    pmids_doi[record["ArticleIds"]["pubmed"][0]] = ""
        except RuntimeError as e:
            print("Warning: Cannot find PMID {} in PUBMED!".format(str(e).split(":")[0]))
            invalid_pmid += 1
        handle.close()
    print("Warnin: {} PMIDs cannot be found in PUMBED, please check!".format(invalid_pmid))
    return pmids_doi


def output_pmid(pmids_doi, output_file):
    """
    Output all pmid-doi pairs to a CSV.
    Each pmid takes one row.
    :param pmids_doi: pmid-doi dictionary;
    :param output_file: path of output file;
    :return: None.
    """
    with open(output_file, "w", newline="") as output:
        print("Total # of PMIDs: {}".format(len(pmids_doi)))
        count_pmids_dois = set(pmids_doi.values())
        count_pmids_dois.remove("")  # remove empty DOIs from count
        print("Total # of PMID-DOI pairs: {}".format(len(count_pmids_dois)))
        csv_writer = csv.writer(output)
        csv_writer.writerow(["PMID", "DOI"])
        [csv_writer.writerow([pmid, doi]) for pmid, doi in pmids_doi.items()]


def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    email = sys.argv[3]
    pmids = get_pmid(input_file)
    pmids_doi = get_doi_api(pmids, email)
    output_pmid(pmids_doi, output_file)
    # print("Total count of pmid: {}".format(len(pmids)))


if __name__ == "__main__":
    main()
