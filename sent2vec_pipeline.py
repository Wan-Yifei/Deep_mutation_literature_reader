import os
import pickle
import sent2vec
from hmac_signature import *
from fetch_pmid_title_abstract import pickle_output


def load_model(model_path):
    """
    Load pre-trained sent2vec model.

    :param model_path: pathway of pre-trained model;
    :return: sen2tev model object.
    """
    model = sent2vec.Sent2vecModel()
    try:
        print("Loading pre-trained model...")
        model.load_model(model_path)
    except Exception as e:
        print(e)
    print('model successfully loaded')
    return model


def embed_sentences(model, sentences_list):
    """
    Embed sentences in a list to sentence_vectors.
    :param model: sent2vec model;
    :param sentences_list: list includes tokenized sentences;
    :return: vectors of sentences in a list.
    """
    sentence_vector = model.embed_sentence(sentences_list)
    return sentence_vector


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


def main():
    model_path = os.environ["MODEL_LOCATION"]
    input_path = os.environ["RAW_CONTENT_PATH"]
    try:
        hmac_key = os.environ["HMAC_KEY"]
    except Exception as e:
        print(e)
        print("Did not receive key of Hash message authentication.")
        hmac_key = ""
    raw_content = input_preprocess(input_path, hmac_key)
    sentences_vetors = {}
    model = load_model(model_path)
    for pmid, sentences in raw_content.values():
        print("\nEmbedding PMID: {}...".format(sentences))
        sentences_vetors = embed_sentences(model, sentences)  # sentences is a list
        sentences_vetors[pmid] = sentences_vetors
    output_file = "embedded_pmid_text.pickle"
    output_folder = "."
    pickle_output(sentences_vetors, output_file, output_folder, key="", use_hmac=False)


if __name__ == '__main__':
    main()
