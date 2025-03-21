from main import find_tandem_repeats_idx, find_repeating_pattern


def test_find_repeating_pattern():
    assert find_repeating_pattern("ATATAT", 2) == "AT"
    assert find_repeating_pattern("ATCG", 2) == "ATCG"
    assert find_repeating_pattern("GCGCGCGC", 2) == "GC"


def test_find_tandem_repeats_idx():

    assert find_tandem_repeats_idx("ATATATAT", min_length=2, min_repeats=3) == [(0, 8)]

    assert find_tandem_repeats_idx("ATCGATCG", min_length=2, min_repeats=3) == []

    assert find_tandem_repeats_idx(
        "AGTTGTTATCGATCTACGGAAGAATT", min_length=2, min_repeats=1
    ) == [(1, 7), (18, 24)]
