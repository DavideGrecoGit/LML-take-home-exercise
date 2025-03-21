import pytest
import json
from main import read_dna_sequences, DNA_SEQUENCE_KEY


def test_valid_file(tmpdir):

    valid_data = {DNA_SEQUENCE_KEY: ["ATCG", "GCTA"]}

    json_file = tmpdir.join("valid_sequences.json")
    json_file.write(json.dumps(valid_data))

    result = read_dna_sequences(str(json_file), DNA_SEQUENCE_KEY)
    assert result == ["ATCG", "GCTA"]


def test_invalid_key(tmpdir):

    valid_data = {DNA_SEQUENCE_KEY: ["ATCG", "GCTA"]}

    json_file = tmpdir.join("valid_sequences.json")
    json_file.write(json.dumps(valid_data))

    with pytest.raises(Exception):
        read_dna_sequences(str(json_file), "other_key")


def test_missing_key(tmpdir):

    invalid_data = {"other_key": ["ATCG", "GCTA"]}

    json_file = tmpdir.join("invalid_sequences.json")
    json_file.write(json.dumps(invalid_data))

    with pytest.raises(Exception):
        read_dna_sequences(str(json_file), DNA_SEQUENCE_KEY)


def test_invalid_file_path():

    invalid_path = "invalid_path.txt"
    with pytest.raises(Exception):
        read_dna_sequences(invalid_path, DNA_SEQUENCE_KEY)


def test_non_json_file(tmpdir):

    text_file = tmpdir.join("non_json_file.txt")

    with pytest.raises(Exception):
        read_dna_sequences(text_file, DNA_SEQUENCE_KEY)
