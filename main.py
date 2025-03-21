import json
import argparse
from pathlib import Path

DNA_SEQUENCE_KEY = "sequences"


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    file_path = parser.add_argument("-file", "--file_path", type=str)
    args = parser.parse_args()

    dna_sequences = read_dna_sequences(args.file_path, DNA_SEQUENCE_KEY)

    print(dna_sequences)
