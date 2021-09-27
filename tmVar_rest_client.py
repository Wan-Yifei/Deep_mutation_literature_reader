import time
import ratelimit
import urllib
from urllib import request
from bs4 import BeautifulSoup

per_second_req_limits = 5  # request limits per period of API
req_time_period = 1  # time period of API requests in seconds


@ratelimit.sleep_and_retry  # Notice: sleep_and_retry decorator must be used before limits decorator.
@ratelimit.limits(calls=per_second_req_limits, period=req_time_period)
def tmvar_rest_api(pmid, concepts, return_format="biocxml", attempts_limit=3, sleep_between_tries=60):
    """
    Call tmVar REST API to recognize variants from abstract of specified PMID.

    :param pmid: string, PMID of the literature; :param concepts: support five kinds of bioconcepts, i.e., Gene,
    Disease, Chemical, Species, Mutation. When 'BioConcept' is used, all five are included.
    :param concepts: request response content, one of Gene, Disease, Chemical, Species, Mutation;
    :param return_format: string, one of PubTator (tab-delimited text file), BioC (xml), and JSON;
    :param attempts_limit: int, How many times failed requests will be automatically retried on error, default by 3;
    :param sleep_between_tries: int, the delay, in seconds, before retrying a request on error;
    :return: bs4.BeautifulSoup object.
    """
    tmvar_root = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/"
    url = tmvar_root + "{}?pmids={}&concepts={}".format(return_format, pmid, concepts)  # url of tmvar with query
    soup = None
    print("Try to call tmVar API of PMID: {}".format(pmid))
    print("From URL: {}".format(url))
    attempts = 0
    while attempts < attempts_limit:
        try:
            attempts += 1
            print("Call tmVar attempt {}".format(attempts))
            response = urllib.request.urlopen(url)
            response = response.read()
            soup = BeautifulSoup(response, features="html.parser")
            break
        except Exception as e:
            print("Error: {}".format(e))
            print("Query tmVar for PMID: {} failed!!".format(pmid))
            print("Sleep {} seconds before retry!\n".format(sleep_between_tries))
            if attempts < attempts_limit:
                time.sleep(sleep_between_tries)
    return soup


def retrieve_variant_entity(soup):
    """
    Retrieve variants list from BeautifulSoup object;
    :param soup: BeautifulSoup Object;
    :return: entities, list of variants.
    """
    anno_list = soup.find_all("annotation")
    entities = []
    for anno in anno_list:
        entities.append(anno.find("text").text)
    entities = list(set(entities))
    return entities


def test_mutation(request_time=10):
    """
    Integrated test for mutation queries.

    :param request_time: count of query times;
    :return: None.
    """
    for i in range(request_time):
        soup = tmvar_rest_api(pmid="24451227", concepts="mutation")
        entities = retrieve_variant_entity(soup)
        print("{}\n".format(entities))


if __name__ == "__main__":
    test_mutation()
