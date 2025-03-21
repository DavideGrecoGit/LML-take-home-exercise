import json
import argparse
from pathlib import Path
from collections import Counter
import re

FILE_PATH = "./dna_sequences.json"
DNA_SEQUENCE_KEY = "sequences"
VALID_CHARS = {"T", "G", "C", "A"}


def read_dna_sequences(file_path: str, dna_sequence_key: str) -> list:
    """Reads DNA sequences from a JSON file.

    Args:
        file_path (str): Path to the JSON file.
        dna_sequence_key (str): The key in the JSON data that corresponds to the list of DNA sequences.

    Raises:
        Exception: If the file path is not valid.
        Exception: If the expected key is not found in the JSON data or if the value associated with the key is not a list.

    Returns:
        list: A list of DNA sequences if found and valid.
    """

    if Path(file_path).is_file() and Path(file_path).suffix == ".json":

        file = open(file_path)
        dna_sequences = json.load(file)

        if dna_sequence_key in dna_sequences:

            sequences = dna_sequences[dna_sequence_key]

            if isinstance(sequences, list):
                return sequences

        raise Exception(f"Could not find suitable DNA sequences.")

    raise Exception("Given file path is not valid.")


def print_invalid_chars(dna_sequences: list[str], valid_chars: set[str]) -> None:
    """For each sequence, prints all the invalid characters found and their positions.
    Additionally, print the count of each invalid character found.

    Args:
        dna_sequences (list[str]): list of DNA sequences.
        valid_chars (list[str]): set of accepted nucleotide bases.
    """

    invalid_chars = []

    for id, sequence in enumerate(dna_sequences):

        if isinstance(sequence, str):

            for j in range(len(sequence)):

                if sequence[j] not in valid_chars:
                    print(f"Sequence {id} - position {j}: {sequence[j]}")
                    invalid_chars.append(sequence[j])
        else:
            print(f"Sequence {id} is not a string")

    invalid_chars = Counter(invalid_chars)

    print("\nInvalid characters: ")
    if len(invalid_chars) > 0:
        for char, count in invalid_chars.items():
            print(f"Char {char} - Count {count}")
    else:
        print("None found")


def remove_invalid_chars(dna_sequences: list[str], valid_chars: set[str]) -> list[str]:
    """For each sequence, remove all the invalid characters found.

    Args:
        dna_sequences (list[str]): list of DNA sequences.
        valid_chars (list[str]): set of accepted nucleotide bases.
    Returns:
        list[str]: A list of DNA sequences without any invalid characters.
    """

    valid_pattern = "".join(valid_chars)

    if valid_pattern == "":
        return ["" for i in dna_sequences]

    for i in range(len(dna_sequences)):

        dna_sequences[i] = re.sub(f"[^{valid_pattern}]", "", dna_sequences[i])

    return dna_sequences


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    file_path = parser.add_argument("-file", "--file_path", type=str, default=FILE_PATH)
    seq_key = parser.add_argument(
        "-key", "--seq_key", type=str, default=DNA_SEQUENCE_KEY
    )
    args = parser.parse_args()

    dna_sequences = read_dna_sequences(args.file_path, args.seq_key)

    print("\n~ Checking for invalid characters in sequence\n")
    print_invalid_chars(dna_sequences, VALID_CHARS)

    print("\n~ Removing invalid characters")
    dna_sequences = remove_invalid_chars(dna_sequences, VALID_CHARS)
    print_invalid_chars(dna_sequences, VALID_CHARS)
