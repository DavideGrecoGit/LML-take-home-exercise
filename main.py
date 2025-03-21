import json
import argparse
from pathlib import Path
from collections import Counter
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

FILE_PATH = "./dna_sequences.json"
DNA_SEQUENCE_KEY = "sequences"
VALID_CHARS = {"T", "G", "C", "A"}

PLOTS_PATH = "./plots"


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
        dna_sequences (list[str]): List of DNA sequences.
        valid_chars (list[str]): Set of accepted nucleotide bases.
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
        dna_sequences (list[str]): List of DNA sequences.
        valid_chars (list[str]): Set of accepted nucleotide bases.
    Returns:
        list[str]: A list of DNA sequences without any invalid characters.
    """

    valid_pattern = "".join(valid_chars)

    if valid_pattern == "":
        return ["" for i in dna_sequences]

    for i in range(len(dna_sequences)):

        dna_sequences[i] = re.sub(f"[^{valid_pattern}]", "", dna_sequences[i])

    return dna_sequences


def compute_gc_content(dna_sequences: list[str]) -> pd.DataFrame:
    """Compute the overall GC-content percentage of each sequence in a given list of DNA sequence.
    The GC-content percentage is calculated as: (number of G and C) / (length of the sequence) * 100

    Args:
        dna_sequences (list[str]): List of DNA sequence, should contain only A, T, C, G characters.

    Returns:
        pd.DataFrame: A DataFrame containing the GC-content of each sequence.
    """

    df_gc_content = []

    for id, sequence in enumerate(dna_sequences):

        if len(sequence) == 0:
            gc_content = 0
        else:
            g_count = sequence.count("G")
            c_count = sequence.count("C")

            gc_content = (g_count + c_count) / len(sequence) * 100

        df_gc_content.append({"gc_content": gc_content})

    df_gc_content = pd.DataFrame(df_gc_content)
    df_gc_content.rename_axis("sequence_id", inplace=True)

    return df_gc_content


def plot_gc_content(
    df_gc_content: pd.DataFrame,
    plot_name_path: str | os.PathLike = None,  # type: ignore
):
    """Compute the overall GC-content of a given list of DNA sequences.
    If a path is given, plot and saves an histogram of the computed gc-content.

    Args:
        df_gc_content (pd.DataFrame): A 1-D DataFrame containing the GC-content of each sequence.
        plot_name_path (str | os.PathLike, optional): Path to save the plot. Defaults to None.

    """

    sns.histplot(data=df_gc_content, kde=True)
    plt.title("Distribution of overall GC Content across DNA Sequences")
    plt.xlabel("GC-content (%)")
    plt.savefig(plot_name_path)
    plt.close()


def count_kmers(sequence: str, k: int = 2) -> dict[str, int]:
    """Counts the k-mers of lenght k contained in the given DNA sequence. Overlapping k-mers are included.

    A k-mer is a sequence of k nucleotides in a DNA sequence.

    Args:
        sequence (str): DNA sequence, should contain only A, T, C, G characters.
        k (int, optional): the k-mer length. Defaults to 2.

    Returns:
        dict[str, int]: dictionary having the found k-mers (keys) and their relative count (values).
    """

    kmers = {}

    if not sequence or k <= 0:
        return kmers

    for i in range(len(sequence) - k + 1):

        key = sequence[i : i + k]
        kmers[key] = kmers.get(key, 0) + 1

    return kmers


def compute_dinucleotide_freq(dna_sequences: list[str]) -> pd.DataFrame:
    """Compute the dinucleotide frequency of a given list of DNA sequences.
    For each found dinucleotide in a sequence, its frequency is computed as: (occurrence of dinucleotide) / (total number of dinucleotides)

    Args:
        dna_sequences (list[str]): List of DNA sequences, should contain only A, T, C, G characters.

    Returns:
        pd.DataFrame: A Pandas DataFrame containing the dinucleotide frequencies for each sequence.
    """

    dinucleotide_frequencies = []

    for sequence in dna_sequences:

        dinucleotide_counts = count_kmers(sequence, 2)

        total_dinucleotides = sum(dinucleotide_counts.values())
        dinucleotide_frequencies.append(
            {
                key: value / total_dinucleotides
                for key, value in dinucleotide_counts.items()
            }
        )

    df_dinucleotide_freq = pd.DataFrame(dinucleotide_frequencies)
    df_dinucleotide_freq.fillna(0, inplace=True)
    df_dinucleotide_freq.rename_axis("sequence_id", inplace=True)

    return df_dinucleotide_freq


def plot_kmer_heatmap(
    df_kmer: pd.DataFrame, plot_title: str, plot_name_path: str | os.PathLike
) -> None:
    """Plot a heatmap of k-mer frequencies.

    Args:
        df_kmer (pd.DataFrame): A DataFrame having as columns all the found kmers and as values their count or frequency in a specific sequence.
        plot_title (str): Title for the plot.
        plot_name_path (str | os.PathLike): Path to save the plot.
    """

    sns.heatmap(df_kmer)
    plt.title(plot_title)
    plt.xlabel("k-mers")
    plt.savefig(plot_name_path)
    plt.close()


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

    print("\n~ a. Calculate and report basic sequence statistic")

    print("Compute and plot overall GC-content")
    df_gc_content = compute_gc_content(dna_sequences)
    plot_gc_content(
        df_gc_content,
        os.path.join(PLOTS_PATH, "overall_gc_content.jpeg"),
    )

    print("Compute dinucleotide frequency")
    df_dinucleotide_freq = compute_dinucleotide_freq(dna_sequences)

    plot_kmer_heatmap(
        df_dinucleotide_freq,
        "Dinucleotide Frequency Heatmap",
        os.path.join(PLOTS_PATH, "dinucleotide_freq.jpeg"),
    )
