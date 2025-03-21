from main import get_complement, find_palindromes_idx


def test_get_complement():
    assert get_complement("A") == "T"
    assert get_complement("T") == "A"
    assert get_complement("C") == "G"
    assert get_complement("G") == "C"
    assert get_complement("E") == ""
    assert get_complement("") == ""


def test_find_palindromes_idx_input_with_valid_palidromes():
    assert find_palindromes_idx("ATTTGCAAAT", min_length=6) == [(0, 10)]
    assert find_palindromes_idx("GGGATTGCAATAAA", min_length=3) == [(3, 11)]
    assert find_palindromes_idx("GGGATTGCAATAAAGCCCAATTGGGCT", min_length=3) == [
        (3, 11),
        (13, 27),
    ]


def test_find_palindromes_idx_input_without_valid_palidromes():
    assert find_palindromes_idx("GGGGGGG", min_length=6) == []
    assert find_palindromes_idx("", min_length=6) == []


def test_find_palindromes_idx_min_lenght_out_of_range():
    assert find_palindromes_idx("ATTTGCAAAT", min_length=30) == []
    assert find_palindromes_idx("ATTTGCAAAT", min_length=-1) == []


def test_find_palindromes_idx_not_enough_bases():
    assert find_palindromes_idx("ATTTTAAAAT", min_length=2) == []
