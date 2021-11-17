import sys
import csv
from pickle_toolbox import *


def retrieve_pmid_title(input_file):
    with open(input_file, "r", encoding="utf-8") as csvfile:
        raw = csv.reader(csvfile)
        header = next(raw)
        header = [colname.lower() for colname in header]
        pmid_ind = header.index("pmid")
        title_ind = header.index("title")
        pmids_title = {row[pmid_ind]: row[title_ind] for row in raw}
    return pmids_title


def retrieve_pmid_source(input_file):
    with open(input_file, "r", encoding="utf-8") as csvfile:
        raw = csv.reader(csvfile)
        header = next(raw)
        header = [colname.lower() for colname in header]
        pmid_ind = header.index("pmid")
        source_ind = header.index("source")
        pmid_source = {}
        for row in raw:
            pmid_source[row[pmid_ind]] = list(set(pmid_source.get(row[pmid_ind], []) + [row[source_ind]]))
    return pmid_source


def main():
    input_file = sys.argv[1]
    output_pmid_title = sys.argv[2]
    output_pmid_source = sys.argv[3]
    try:
        output_hmac_key = sys.argv[4]
    except IndexError:
        print("Notice: No HMAC key provided!")
        output_hmac_key = None
    pmid_title = retrieve_pmid_title(input_file)
    pmid_source = retrieve_pmid_source(input_file)
    if output_hmac_key:
        pickle_output(pmid_title, output_pmid_title, output_folder=".", key=output_hmac_key, use_hmac=True)
        pickle_output(pmid_source, output_pmid_source, output_folder=".", key=output_hmac_key, use_hmac=True)
    else:
        pickle_output(pmid_title, output_file=output_pmid_title, output_folder=".")
        pickle_output(pmid_source, output_file=output_pmid_source, output_folder=".")


if __name__ == '__main__':
    main()
