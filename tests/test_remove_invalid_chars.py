from main import remove_invalid_chars, VALID_CHARS


def test_valid_sequences():

    dna_sequences = ["ATCG", "GCTA", "CCCCAAAA"]

    result = remove_invalid_chars(dna_sequences, VALID_CHARS)
    assert result == ["ATCG", "GCTA", "CCCCAAAA"]


def test_invalid_characters():

    dna_sequences = ["ATCG", "GCTZ", "XFH", "A1T2C3G"]

    result = remove_invalid_chars(dna_sequences, VALID_CHARS)
    print(result)
    assert result == ["ATCG", "GCT", "", "ATCG"]


def test_empty_input():

    dna_sequences = []

    result = remove_invalid_chars(dna_sequences, VALID_CHARS)
    assert result == []


def test_no_valid_characters():
    VALID_CHARS = set()
    dna_sequences = ["ATCG", "GCTA"]

    result = remove_invalid_chars(dna_sequences, VALID_CHARS)
    assert result == ["", ""]
