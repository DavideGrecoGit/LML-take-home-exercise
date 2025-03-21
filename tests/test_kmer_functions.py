import pandas as pd
from main import count_kmers, compute_dinucleotide_freq, compute_top_kmers_count


def test_count_kmers_valid_input():
    assert count_kmers("AAAA", 1) == {"A": 4}
    assert count_kmers("AATT", 2) == {"AA": 1, "AT": 1, "TT": 1}
    assert count_kmers("ATCGATCG", 2) == {"AT": 2, "TC": 2, "CG": 2, "GA": 1}
    assert count_kmers("GCAT", 3) == {"GCA": 1, "CAT": 1}
    assert count_kmers("GCGC", 4) == {"GCGC": 1}


def test_count_kmers_empty_input():
    assert count_kmers("", 10) == {}
    assert count_kmers("", 0) == {}


def test_count_kmers_incorrect_k():
    assert count_kmers("A", 75) == {}
    assert count_kmers("ATGC", 0) == {}
    assert count_kmers("ATATGCG", -1) == {}


def test_compute_dinucleotide_freq_valid_sequences():
    dna_sequences = ["AATT", "GCGC"]

    expected_output = pd.DataFrame(
        [
            {"AA": 1 / 3, "AT": 1 / 3, "TT": 1 / 3},
            {"GC": 2 / 3, "CG": 1 / 3},
        ]
    ).fillna(0)
    expected_output.rename_axis("sequence_id", inplace=True)

    result = compute_dinucleotide_freq(dna_sequences)
    pd.testing.assert_frame_equal(result, expected_output)


def test_compute_dinucleotide_freq_empty_input():
    result_empty = compute_dinucleotide_freq([])
    expected_empty = pd.DataFrame().fillna(0)
    expected_empty.rename_axis("sequence_id", inplace=True)

    pd.testing.assert_frame_equal(result_empty, expected_empty)


def test_top_kmers_count_valid_input():
    dna_sequences = ["ATATATCGATG", "GCGC", "AAATTTCGGG"]

    expected_output = pd.DataFrame([{"AT": 4}, {"GC": 2}, {"AA": 2}]).fillna(0)
    expected_output.rename_axis("sequence_id", inplace=True)

    result = compute_top_kmers_count(dna_sequences, k=2, top_n=1)
    pd.testing.assert_frame_equal(result, expected_output)


def test_top_kmers_count_n_higher_than_possible_kmers():

    dna_sequences = ["ATAT", "GCAA", "AA"]

    expected_output = pd.DataFrame(
        [{"AT": 2, "TA": 1}, {"GC": 1, "CA": 1, "AA": 1}, {"AA": 1}]
    ).fillna(0)
    expected_output.rename_axis("sequence_id", inplace=True)

    result = compute_top_kmers_count(dna_sequences, k=2, top_n=100)
    pd.testing.assert_frame_equal(result, expected_output)


def test_top_kmers_count_negative_n():
    dna_sequences = ["ATAT", "GCAA", "AA"]

    result = compute_top_kmers_count(dna_sequences, k=2, top_n=-1)
    expected_output = pd.DataFrame([])
    expected_output.rename_axis("sequence_id", inplace=True)
    pd.testing.assert_frame_equal(result, expected_output)


def test_top_kmers_count_empty_sequences():
    dna_sequences = ["", "", ""]

    expected_output = pd.DataFrame([{}, {}, {}])
    expected_output.rename_axis("sequence_id", inplace=True)

    result = compute_top_kmers_count(dna_sequences, k=2, top_n=1)

    print(result)

    pd.testing.assert_frame_equal(result, expected_output)


def test_top_kmers_count_empty_input():
    dna_sequences = []

    expected_output = pd.DataFrame([])
    expected_output.rename_axis("sequence_id", inplace=True)

    result = compute_top_kmers_count(dna_sequences, k=2, top_n=1)

    print(result)

    pd.testing.assert_frame_equal(result, expected_output)
