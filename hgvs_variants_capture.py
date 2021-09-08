import re
import sys

default_var_pattern = {
    "dna": {
        "substitution_change": r"[cCgG]\.(\d+|\*\d+|-\d+)([+-]\d+)?([GCTAgcta])?>([GCTAgcta])",
        "substitution_no_change": r"[cCgG]\.(\d+|\*\d+|-\d+)([+-]\d+)?([GCTAgcta])?=",
        "deletion": r"[cCgG]\.(\d+|\*\d+|-\d+)([+-]\d+)?(?:_(\d+|\*\d+|-\d+)([+-]\d+)?)?del(?!ins)",
        "insertion_nucleotides": r"[cCgG]\.(\d+|\*\d+|-\d+)([+-]\d+)?_(\d+|\*\d+|-\d+)([+-]\d+)?ins([GCTAgcta]+)",
        "insertion_copy": r"[cCgG]\.(\d+|\*\d+|-\d+)([+-]\d+)?_(\d+|\*\d+|-\d+)([+-]\d+)?ins(\d+|\*\d+|-\d+)(["
                          r"+-]\d+)?_(\d+|\*\d+|-\d+)([+-]\d+)?",
        "indel": r"[cCgG]\.(\d+|\*\d+|-\d+)([+-]\d+)?(?:_(\d+|\*\d+|-\d+)([+-]\d+)?)?([GCTAgcta]+)?delins(["
                 r"GCTAgcta]+)",
        "inversion": r"[cCgG]\.(\d+|\*\d+|-\d+)([+-]\d+)?_(\d+|\*\d+|-\d+)([+-]\d+)?inv",
        "duplication": r"[cCgG]\.(\d+|\*\d+|-\d+)([+-]\d+)?(?:_(\d+|\*\d+|-\d+)([+-]\d+)?)?dup([GCTAgcta]+)?"
    }
}


def map_temp_encoding(key, key_2_placeholder, placeholder_prefix, notice):
    """
    Assign a placeholder (temporary encoding) as a value for input key and then update the provided dictionary.

    :param key: string, a key of the dictionary key_2_placeholder;
    :param key_2_placeholder: dictionary, stores key and corresponding temporary placeholder;
    :param placeholder_prefix: string, prefix of temporary placeholders;
    :param notice: bool, if True, print out the dictionary update process for each key;
    :return: None.
    """
    index = len(key_2_placeholder)
    value = placeholder_prefix + str(index)
    if key not in key_2_placeholder:
        key_2_placeholder[key] = value
    if notice:
        print("key {} has been mapped to temporary encoding {}!".format(key, key_2_placeholder[key]))


def recognize_variant_entity(variant_type, text, pattern=default_var_pattern):
    """
    Recognize variant entities in the text by matching regular expression pattern.

    :param variant_type: string, should be one of "dna", "rna", and "protein";
    :param text: string, source text;
    :param pattern: dictionary, regex patterns for each variant type;
    :return: recognized_entities, dictionary: key -> variant type, value -> generator of regex match results.
    """
    recognized_entities = {}  # store all match words
    for vartype, regex in pattern[variant_type].items():
        recognized_entities[vartype] = re.finditer(regex, text)
    return recognized_entities


def build_var_temp_encoding(recognized_entities, placeholder_prefix="temp_placeholder_", notice=False):
    """
    Map variant to placeholder (temporary encoding).

    :param recognized_entities: dictionary, key -> variant type, value -> generator of regex match results;
    :param placeholder_prefix: string, prefix of placeholder, by default "temp_placeholder_";
    :param notice: bool, if True, print out the dictionary update process for each key;
    :return: key_2_placeholder, dictionary: key -> variant, value -> placeholder.
    """
    key_2_placeholder = {}  # dict: variant -> placeholder
    for variants in recognized_entities.values():
        [map_temp_encoding(var.group(0), key_2_placeholder, placeholder_prefix, notice) for var in variants]
    return key_2_placeholder


def replace_var_2_temp_encoding(key_2_placeholder, text, notice=False):
    """
    Replace variant entities in provided text with placeholders.

    :param key_2_placeholder: dictionary: key -> variant, value -> placeholder;
    :param text: string, source text;
    :param notice: bool, if True, print out the replace process for each variant;
    :return: text, string: replaced variants with placeholders.
    """
    for var, placeholder in key_2_placeholder.items():
        if notice:
            print("========================================")
            print("Variant: {}".format(var))
            print("Placeholder: {}\n".format(placeholder))
            print("Original text:\n{}\n".format(text))
        text = text.replace(var, placeholder)
        if notice:
            print("Updated text:\n{}\n".format(text))
    return text


def replace_temp_encoding_2_var(key_2_placeholder, text, notice=False):
    """
    Replace placeholders in provided text with variant entities.

    :param key_2_placeholder: dictionary: key -> variant, value -> placeholder;
    :param text: string, source text;
    :param notice: bool, if True, print out the replace process for each variant;
    :return: text, string: replaced placeholders with variants.
    """
    placeholder_2_key = {placeholder: var for var, placeholder in key_2_placeholder.items()}  # exchange key with value
    for placeholder, var in placeholder_2_key.items():
        if notice:
            print("========================================")
            print("Placeholder: {}\n".format(placeholder))
            print("Variant: {}".format(var))
            print("Original text:\n{}\n".format(text))
        text = text.replace(placeholder, var)
        if notice:
            print("Updated text:\n{}\n".format(text))
    return text


def test(build_notice, key2pla_notice, pla2key_notice):
    """
    Integrated test.

    :param build_notice: bool, if True, print out the dictionary update process for each key;
    :param key2pla_notice: bool, if True, print out the replace process for each variant;
    :param pla2key_notice: bool, if True, print out the replace process for each variant;
    :return: None
    """
    text = "Try regex on DNA level substitution change variants: c.123+45A>G, c.124-56C>T, c.*32G>A. \nAnd " \
           "substitution no change variant c.-123-64G=. \nAlso deletion c.4072-1234_5155-246del. " \
           "\nMoreover, insertion g.32867861_32867862insT and c.849_850ins858_895. \n" \
           "Further, indel c.6775_6777delinsC. \nThen, inversion g.32361330_32361333inv. \n" \
           "Don't forget inversion g.32361330_32361333inv and c.*77-10_*77-1inv. \n" \
           "Beside duplication, c.260_264+48dup."
    recognized_entities = recognize_variant_entity("dna", text)
    key_2_placeholder = build_var_temp_encoding(recognized_entities, notice=build_notice)
    new_text = replace_var_2_temp_encoding(key_2_placeholder, text, notice=key2pla_notice)
    replace_temp_encoding_2_var(key_2_placeholder, new_text, notice=pla2key_notice)


if __name__ == "__main__":
    try:
        print_build = sys.argv[1]
        print_key2pla = sys.argv[2]
        print_pla2key = sys.argv[3]
    except IndexError:
        print_build = False
        print_key2pla = False
        print_pla2key = True
    test(print_build, print_key2pla, print_pla2key)
