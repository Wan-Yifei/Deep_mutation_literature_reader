import sys
import csv
import math
import numpy as np
import matplotlib.pyplot as plt

# binary label of ACMG classification
binary_label_ref = {"neutral": ["BENIGN", "LIKELY BENIGN", "VOUS"], "pathogenic": ["LIKELY PATHOGENIC", "PATHOGENIC"]}


def get_label(source):
    with open(source, "r", encoding='utf-8') as raw:
        csv_raw = csv.DictReader(raw)
        pmid_acmg_label = {row["pmids"]: row["classification"] for row in csv_raw if row["pmids"]}
    return pmid_acmg_label


def dualize_label(raw_label, binary_label=binary_label_ref):
    binary_label_mapping = {}
    for binary_label, acmg_label in binary_label.items():
        for acmg in acmg_label:
            binary_label_mapping[acmg] = binary_label
    pmids_bin_tmp = {pmids: binary_label_mapping[acmg_label] for pmids, acmg_label in raw_label.items()}
    pmid_bin_label = {}
    for pmids, label in pmids_bin_tmp.items():
        pmids = set(pmids.split("|"))
        for pmid in pmids:
            if pmid not in pmid_bin_label.keys():
                pmid_bin_label[pmid] = label
            elif pmid_bin_label[pmid] != label:
                print(pmid)
                print(label)
                raise Exception("Doesn't match!!")
            else:
                pmid_bin_label[pmid] = label
    return pmid_bin_label


def get_label_list(bin_label):
    pmid_list = [pmid for pmid in bin_label.keys()]
    pmid_list = np.array(pmid_list)
    np.random.shuffle(pmid_list)
    label_list = [1 if bin_label[pmid] == "pathogenic" else 0 for pmid in pmid_list]
    label_list = np.array(label_list)
    return label_list, pmid_list

def main():
    source_path = sys.argv[1]
    pmid_acmg_label = get_label(source_path)
    pmid_bin_label = dualize_label(pmid_acmg_label)
    label_list, pmid_list = get_label_list(pmid_bin_label)
    print(label_list.shape)
    print(pmid_list.shape)
    print(pmid_list)
    split_index = math.floor(label_list.shape[0] * 0.98)
    #print("# of 1: {}".format(label_list.count(1)))
    #print("# of 0: {}".format(label_list.count(0)))
    #plt.hist(label_list, bins=[-.5, .5, 1.5])
    #plt.xticks((0, 1))
    #plt.show()
    #plt.close()


if __name__ == "__main__":
    main()
