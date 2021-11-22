import os
import sys
import pickle
import fasttext 
from hmac_signature import *
from fetch_pmid_title_abstract import pickle_output


def load_model(model_path):
    """
    Load pre-trained sent2vec model.

    :param model_path: pathway of pre-trained model;
    :return: sen2tev model object.
    """
    try:
        print("Loading pre-trained model...")
        model = fasttext.load_model(model_path)
    except Exception as e:
        print(e)
        raise Exception("Model loading failed!")
    print('model successfully loaded')
    return model


def embed_words(model, word_list):
    """
    Embed words in a list to word_vectors.
    :param model: fasttext model;
    :param word_list: list includes tokenized words;
    :return: vectors of words in a list.
    """
    word_vectors = [model[word] for word in word_list]
    return word_vectors


def input_preprocess(input_content_path, key=None):
    """
    Read raw sentences from a pickle file.

    :param input_content_path: path of input pickle file;
    :param key: HMAC key;
    :return: content from pickle.
    """
    with open(input_content_path, "rb") as text:
        if key:
            text_signature_old = text.readline().decode("utf-8")
            text_content = b"".join(text.readlines())
            if verify_hmac_signature(text_content, key, text_signature_old):
                raw_content = pickle.loads(text_content)  # A dict
                print("Raw content is ready!")
        else:
            raw_content = pickle.load(text)
    return raw_content


def word_tokenize(sentences_list):
    """
    Tokenize raw sentences for each pmid.

    :param sentences_list: sentences in a list or just a string;
    :return: tokenized words in a list.
    """
    if isinstance(sentences_list, list):
        new_list = " ".join(sentences_list)
        word_token = new_list.split(" ")
        word_token = [word for word in word_token if word]  # remove space from list
    else:
        word_token = sentences_list.split(" ")
        word_token = [word for word in word_token if word]  # remove space from list

    return word_token


def main():
    output_file = sys.argv[1]
    model_path = os.environ["MODEL_LOCATION"]
    input_path = os.environ["RAW_CONTENT_PATH"]
    try:
        hmac_key = os.environ["HMAC_KEY"]
    except Exception as e:
        print(e)
        print("Did not receive key of Hash message authentication.")
        hmac_key = ""
    raw_content = input_preprocess(input_path, hmac_key)
    words_vectors_dict = {}
    model = load_model(model_path)
    n = 0
    for pmid, sentences in raw_content.items():
        if isinstance(sentences, list):
            if sentences[0]:
                n += 1
                print("\n{}-1. Tokenize sentences of PMID: {}".format(n, pmid))
                word_token = word_tokenize(sentences)
                print("{}\n".format(word_token))
                print("\n{}-2. Embedding PMID: {}...".format(n, pmid))
                word_vectors = embed_words(model, word_token)  # word_token is a list
                assert len(word_token) == len(word_vectors), "Length doesn't match."
                words_vectors_dict[pmid] = word_vectors
        else:
            if sentences:
                n += 1
                print("\n{}-1. Tokenize sentences of PMID: {}".format(n, pmid))
                word_token = word_tokenize(sentences)
                print("{}\n".format(word_token))
                print("\n{}-2. Embedding PMID: {}...".format(n, pmid))
                word_vectors = embed_words(model, word_token)  # word_token is a list
                assert len(word_token) == len(word_vectors), "Length doesn't match."
                words_vectors_dict[pmid] = word_vectors

    output_folder = "."
    pickle_output(words_vectors_dict, output_file, output_folder, key="", use_hmac=False)


if __name__ == '__main__':
    main()
