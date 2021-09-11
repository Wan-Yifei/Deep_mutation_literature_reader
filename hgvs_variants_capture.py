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
    key = key.strip(" ,.)")  # remove following space, comma, period and parentheses matched by \s.
    if key not in key_2_placeholder:
        key_2_placeholder[key] = value
    if notice:
        print("key {} has been mapped to temporary encoding {}!".format(key, key_2_placeholder[key]))


def escape_character_variant(variant):
    """
    Add backslash to special characters to convert a variant name string into a regex pattern.

    :param variant: string, variant;
    :return: variant_pattern: regex variant pattern.
    """
    special_characters = '''{}[]()^$.|*+?\\'''  # special chars from string.punctuation
    var_characters = []
    for char in variant:
        if char in special_characters:
            var_characters.append(char.replace(char, r"\{}".format(char)))
        else:
            var_characters.append(char)
    variant_pattern = "".join(var_characters) + r".*?(?:\s|$)"
    variant_pattern = variant_pattern.lower()
    return variant_pattern


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


def get_variant_entity_tmvar(text, ner):
    assert (isinstance(ner, list)), "Please check your ner arguments."
    recognized_entities = {}  # store all match words
    n = 0
    # the ner list includes variants entities recognized by other tool
    for variant in ner:
        regex = escape_character_variant(variant)
        vartype = "tmvar_" + str(n)
        recognized_entities[vartype] = re.finditer(regex, text)  # variant entities from tmVar tool
        n += 1
    return recognized_entities


def build_var_temp_encoding(recognized_entities, key_2_placeholder, placeholder_prefix="placeholder", notice=False):
    """
    Map variant to placeholder (temporary encoding).

    :param recognized_entities: dictionary, key -> variant type, value -> generator of regex match results;
    :param key_2_placeholder: dictionary, key -> variant, value -> placeholder, store mapping between key/placeholder;
    :param placeholder_prefix: string, prefix of placeholder, by default "temp_placeholder_";
    :param notice: bool, if True, print out the dictionary update process for each key;
    :return: key_2_placeholder, dictionary: key -> variant, value -> placeholder.
    """
    for key, variants in recognized_entities.items():
        [map_temp_encoding(var.group(0), key_2_placeholder, placeholder_prefix, notice) for var in variants]
    if notice:
        print(key_2_placeholder)
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
            print("Placeholder: {}".format(placeholder))
            print("Variant: {}\n".format(var))
            print("Original text:\n{}\n".format(text))
        pattern = r"{}\b".format(placeholder)
        text = re.sub(pattern, var, text)
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


def test2(build_notice, key2pla_notice, pla2key_notice):
    """
    Integrated test2. Test both regex method and tmVar method.

    :param build_notice: bool, if True, print out the dictionary update process for each key;
    :param key2pla_notice: bool, if True, print out the replace process for each variant;
    :param pla2key_notice: bool, if True, print out the replace process for each variant;
    :return: None
    """
    text = "c.224g>a, p.arg75gln (r75q) presumably leads to an amino-acid change from arginine to glutamine in the " \
           "membrane-spanning domain of the cftr protein. initially reported as a benign sequence variation, " \
           "p.arg75gln was associated with a high risk of pancreatitis, a risk that was strikingly higher " \
           "when p.arg75gln was combined with a spink1 variant. in addition, it was shown that p.arg75gln alters " \
           "bicarbonate but not chloride conductance and that the mutation also induces exon 3 skipping. to " \
           "investigate the role of p.arg75gln in idiopathic chronic pancreatitis (icp), we performed genotyping of " \
           "the cftr gene in 880 patients with icp, 198 patients with idiopathic bronchiectasis (ib), 74 patients " \
           "with classical cystic fibrosis (cf), 48 patients with congenital bilateral absence of the vas deferens (" \
           "cbavd) and 148 healthy controls. p.arg75gln was identified in 3.3% (29/880) of patients with icp, " \
           "3.3% (9/272) patients with a pulmonary disease, 2.1% (1/48) of patients with cbavd and 4.7% (7/148) of " \
           "healthy controls. it was frequently associated with the c.[1210-12t[7];1408a>g] (t7-p.val470) allele and " \
           "this cftr genetic background could not explain the putative pathogenicity of this variant. to assess " \
           "whether cftr and spink1 mutations are co-inherited in pancreatitis, we sequenced spink1 gene exon 3 in " \
           "the 46 patients who were previously identified to be heterozygous for p.arg75gln. two spink1 " \
           "pancreatitis-associated variants, p.asn34ser and p.pro55ser, were found in 6 patients: 4 of 29 (13.8%) " \
           "patients with icp (3 p.asn34ser and 1 p.pro55ser), 1 of 7 (14.3%) healthy controls (p.asn34ser) and 1 of " \
           "9 (11.1%) patients with ib (p.pro55ser). our study does not confirm that the cftr p.arg75gln mutation " \
           "confers a significant risk of pancreatitis when considered individually and with a concurrent spink1 " \
           "mutation, suggesting the role of other genetic and environmental factors. "
    tmvar_ner = ["p.Arg75Gln", "p.Pro55Ser", "c.[1210-12T[7]", "p.Asn34Ser", "R75Q"]
    print("Original text: \n{}".format(text))
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    tmvar_entities = get_variant_entity_tmvar(text, tmvar_ner)
    key_2_placeholder_init = {}
    # Build tmVar based key to placeholder mapping
    key_2_placeholder_tmvar = build_var_temp_encoding(tmvar_entities, key_2_placeholder_init, notice=build_notice)
    tmvar_text = replace_var_2_temp_encoding(key_2_placeholder_tmvar, text, notice=key2pla_notice)  # tmVar replacement
    print("Text after tmVar replacement: \n{}".format(tmvar_text))
    regex_entities = recognize_variant_entity("dna", tmvar_text)  # conventional regex NER
    # Combine regex based key2placeholder with tmVar based key2placeholder
    key_2_placeholder = build_var_temp_encoding(regex_entities, key_2_placeholder_tmvar, notice=build_notice)
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    regex_tmvar_text = replace_var_2_temp_encoding(key_2_placeholder, tmvar_text, notice=key2pla_notice)
    new_text = replace_temp_encoding_2_var(key_2_placeholder, regex_tmvar_text, notice=pla2key_notice)
    print("Text after both tmVar and regex based replacement: \n{}".format(regex_tmvar_text))
    assert (text == new_text), "Error: New text does not match the original text."
    print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
    print("\nNew text matches the original text! Test passed!")


if __name__ == "__main__":
    try:
        print_build = sys.argv[1]
        print_key2pla = sys.argv[2]
        print_pla2key = sys.argv[3]
    except IndexError:
        print_build = False
        print_key2pla = False
        print_pla2key = False
    test2(print_build, print_key2pla, print_pla2key)
