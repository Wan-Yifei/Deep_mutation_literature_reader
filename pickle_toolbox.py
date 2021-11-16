import pickle
from hmac_signature import *


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


def pickle_load(pickle_path, hmac_key=None, use_hmac=False):
    """
    Load pickle file.

    :param pickle_path: string, the path of pickled file.
    :param hmac_key: string, the key used by HASH signature.
    :param use_hmac: bool, True/False.
    :return: un-pickled content.
    """
    with open(pickle_path, "rb") as pickled_file:
        if use_hmac:
            signature = pickled_file.readline().decode("utf-8")
            content = b"".join(pickled_file.readlines())
            if not hmac_key:
                raise Exception("Key value of HASH message authentication MISSED!!")
            if not verify_hmac_signature(content, hmac_key, signature):
                raise Exception("Key value of HASH message authentication is not CORRECT!!")
            unpickle_file = pickle.loads(content)
        else:
            try:
                unpickle_file = pickle.loads(pickled_file)
            except Exception as e:
                print("Does pickle file use HASH message authentication???\n")
                raise e
        print("The content of the pickled file is ready!")
        return unpickle_file
