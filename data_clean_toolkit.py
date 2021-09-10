import re
import string


def remove_url(text):
    """
    Remove url from input text.

    :param text: string, input text;
    :return: string without url.
    """
    return re.sub(r"https?://\S+|www\.\S+", "", text)


def remove_html(text):
    """ Remove html from input text.

    :param text: string, input text;
    :return: string without html.
    """
    html = re.compile(r"<.*?>|&([a-z0-9]+|#[0-9]{1,6}|#x[0-9a-f]{1,6});")
    return re.sub(html, "", text)


def remove_special_characters(text):
    """
    Remove special special characters, including symbols, emojis, and other graphic characters.

    :param text: string, input text;
    :return: string without special characters.
    """
    emoji_pattern = re.compile(
        '['
        u'\U0001F600-\U0001F64F'  # emoticons
        u'\U0001F300-\U0001F5FF'  # symbols & pictographs
        u'\U0001F680-\U0001F6FF'  # transport & map symbols
        u'\U0001F1E0-\U0001F1FF'  # flags (iOS)
        u'\U00002702-\U000027B0'
        u'\U000024C2-\U0001F251'
        u'\U000024C2-\U0001F251'
        u'\u00A9'                 # copyright symbol
        ']+',
        flags=re.UNICODE)
    return emoji_pattern.sub(r'', text)


def remove_punctuation(text ):
    """
    Remove punctuation from input string.

    :param text: string, input text;
    :param variant_safe: bool, if true, don't change variant name. e.g. keep c.188T>C while true else output c188TC;
    :return: string without punctuation.
    """
    return text.translate(str.maketrans('', '', string.punctuation))
