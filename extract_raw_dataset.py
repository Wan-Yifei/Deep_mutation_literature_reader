import sys
import csv


def retrieve_pmid_title(input_file):
    with open(input_file, "r", encoding="utf-8") as csvfile:
        raw = csv.reader(csvfile)
        header = next(raw)
        header = [colname.lower() for colname in header]
        pmid_ind = header.index("pmid")
        title_ind = header.index("title")
        pmids_title = {row[pmid_ind]: row[title_ind] for row in raw}
    return pmids_title


def retrieve_pmid_label(input_file):
    with open(input_file, "r", encoding="utf-8") as csvfile:
        raw = csv.reader(csvfile)
        header = next(raw)
        header = [colname.lower() for colname in header]
        pmid_ind = header.index("pmid")
        source_ind = header.index("source")
        tmp = {}
        pmid_source = {}
        for row in raw:
            pmid_source[row[pmid_ind]] = list(set(pmid_source.get(row[pmid_ind], []) + [row[source_ind]]))
    return pmid_source


def main():
    input_file = sys.argv[1]
    pmid_title = retrieve_pmid_title(input_file)
    pmid_source = retrieve_pmid_label(input_file)


if __name__ == '__main__':
    main()