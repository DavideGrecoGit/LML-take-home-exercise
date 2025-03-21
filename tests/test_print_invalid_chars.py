import pytest
from main import print_invalid_chars, VALID_CHARS


def test_valid_sequences(capfd):
    valid_chars = {"A", "T", "C", "G"}
    dna_sequences = ["ATCG", "GCTA"]

    print_invalid_chars(dna_sequences, valid_chars)

    captured = capfd.readouterr()
    assert "Invalid characters: " in captured.out
    assert "None found" in captured.out


def test_empty_input(capfd):
    valid_chars = {"A", "T", "C", "G"}
    dna_sequences = []

    print_invalid_chars(dna_sequences, valid_chars)

    captured = capfd.readouterr()
    assert "Invalid characters: " in captured.out
    assert "None found" in captured.out


def test_invalid_characters(capfd):
    valid_chars = {"A", "T", "C", "G"}
    dna_sequences = ["AFFCG", "GCTZ"]

    print_invalid_chars(dna_sequences, valid_chars)

    captured = capfd.readouterr()
    assert "Sequence 0 - position 2: F" in captured.out
    assert "Sequence 0 - position 1: F" in captured.out
    assert "Sequence 1 - position 3: Z" in captured.out
    assert "Invalid characters: " in captured.out
    assert "Char F - Count 2" in captured.out
    assert "Char Z - Count 1" in captured.out


def test_non_string_sequence(capfd):
    valid_chars = {"A", "T", "C", "G"}
    dna_sequences = ["ATCG", 123, "GCTA"]

    print_invalid_chars(dna_sequences, valid_chars)

    captured = capfd.readouterr()
    assert "Sequence 1 is not a string" in captured.out
    assert "Invalid characters: " in captured.out
    assert "None found" in captured.out
