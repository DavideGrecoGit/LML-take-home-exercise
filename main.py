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

K_TOP_MERS = [3, 4, 5]

PALINDROMES_MIN_LENGTH = 20

TANDEM_REPEATS_MIN_LENGTH = 1
TANDEM_REPEATS_MAX_LENGTH = 100
TANDEM_REPEATS_MIN_REPEATS = 5


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
    """Compute the dinucleotide frequency for each sequence in a given list of DNA sequences.
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


def compute_top_kmers_count(
    dna_sequences: list[str], k: int, top_n: int
) -> pd.DataFrame:
    """Compute the counts of the top k-mers for each sequence in a given list of DNA sequences.

    A k-mer is a sequence of k nucleotides in a DNA sequence.

    Args:
        dna_sequences (list[str]): List of DNA sequences, should contain only A, T, C, G characters.
        k (int): the k-mer length.
        top_n (int): The number of top k-mers to return. If negative returns an empty Pandas Dataframe.

    Returns:
        pd.DataFrame: A Pandas DataFrame containing the count of the top n k-mers for each sequence.
    """

    df_top_kmers = []

    if top_n > 0:

        for sequence in dna_sequences:

            kmers_counts = count_kmers(sequence, k)

            top_keys = sorted(kmers_counts, key=kmers_counts.get, reverse=True)[:top_n]  # type: ignore

            df_top_kmers.append({key: kmers_counts[key] for key in top_keys})

    df_top_kmers = pd.DataFrame(df_top_kmers)
    df_top_kmers.fillna(0, inplace=True)
    df_top_kmers.rename_axis("sequence_id", inplace=True)

    return df_top_kmers


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


def get_complement(nucleotide: str) -> str:
    """Get the complementary nucleotide of a given nucleotide base.

    Args:
        nucleotide (str): A single nucleotide (A, T, C, G).

    Returns:
        str: The complementary nucleotide. If no complementary is found, an empty string is returned.
    """

    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    return complement.get(nucleotide, "")


def find_palindromes_idx(sequence: str, min_length: int) -> list[tuple[int, int]]:
    """Find the start and end indices of all palindromes in the given sequence that
    are even in lenght, longer than or equal to the specified minimum length and contain at least three different nucleotide bases.

    A nucleotide sequence is palindromic if it is equal to its reverse complement.

    Args:
        sequence (str): DNA sequence, should contain only A, T, C, G characters.
        min_length (int): Desired minimum length of the palindrome.

    Returns:
        list[tuple[int, int]]: List of all the start and end indices of the found palindromes.
                                If no palindromes are found, an empty list is returned.
    """

    palindromes = []

    if min_length < 3:
        return palindromes

    for i in range(len(sequence) - 1):
        left = i
        right = i + 1

        while (
            left >= 0
            and right < len(sequence)
            and sequence[left] == get_complement(sequence[right])
        ):
            left -= 1
            right += 1

        pal = sequence[left + 1 : right]

        if len(pal) >= min_length and len(set(pal)) >= 3:
            palindromes.append((left + 1, right))

    return palindromes


def collect_all_palindromes(
    dna_sequences: list[str], min_length: int = 20
) -> pd.DataFrame:
    """In each sequence of a given list of DNA sequences, find all palindromes along with:
     - their lenght
     - start and end indices

    Args:
        dna_sequences (list[str]): List of DNA sequences, should contain only A, T, C, G characters.
        min_length (int, Optional): Desired minimum length of the palindrome. Defaults to 20.
    Returns:
        pd.DataFrame: A Pandas DataFrame containing the found palindromes and their information, for each sequence.
    """
    all_palindromes = []

    for id, sequence in enumerate(dna_sequences):

        palindromes_idx = find_palindromes_idx(sequence, min_length=min_length)

        for start, end in palindromes_idx:
            all_palindromes.append(
                {
                    "sequence_id": id,
                    "palindrome": sequence[start:end],
                    "length": end - start,
                    "start": start + 1,  # Adjusting for 1-based index
                    "end": end,
                }
            )

    df_palindromes = pd.DataFrame(all_palindromes)
    df_palindromes.sort_values("length", ascending=False)

    return df_palindromes


def find_repeating_pattern(tandem_repeat: str, min_length: int) -> str:
    """Find the repeating pattern of a given tandem_repeat sequence.
    The repeating pattern must be longer than or equal to the specified minimum length.

    A repeating pattern is the smallest repeating unit that can be repeated to form the given tandem repeat.

    Args:
        tandem_repeat (str): The tandem repeat sequence to decompone.
        min_length (int): Minimum lenght of the repeating pattern.

    Returns:
        str: The found repeating pattern. If no repeating pattern is found, an empty string is returned.
    """

    for i in range(min_length - 1, len(tandem_repeat)):

        base_pattern = tandem_repeat[: i + 1]

        if base_pattern * (len(tandem_repeat) // len(base_pattern)) == tandem_repeat:
            return base_pattern

    return ""


def find_tandem_repeats_idx(
    sequence: str,
    min_length: int = 1,
    max_length: int = -1,
    min_repeats: int = 3,
) -> list[tuple[int, int]]:
    """Find the start and end indices of all tandem repeats in the given sequence that have a length between
    the specified minimum and maximum lengths (inclusive) and meet the minimum number of repetitions required.
    Overlapping tandem repeats are NOT included.

    A tandem repeat is a DNA sequence consisting of a repeating pattern occurring consecutively.

    Args:
        sequence (str): DNA sequence, should contain only A, T, C, G characters.
        min_length (int): The minimum length of the repeating pattern.
        max_length (int, optional): The maximum length of the repeating pattern. Defaults to -1,
                                     which indicates no maximum length.
        min_repeats (int): The minimum number of repeats required to consider it a tandem repeat.

    Returns:
        list[tuple[int, int]]: List of all the start and end indices (0-based) of the found tandem repeats.
                               If no tandem repeats are found, an empty list is returned.

    """

    # (.{min_length,max_lenght}) -> all sequences of any char of length between min and max length
    # \1 -> any matched content has to be repeated
    # {min_repeats,} -> they have to repeat at least "min_repeats" times
    pattern = r"(.{%d,%d})\1{%d,}" % (min_length, max_length, min_repeats)

    if min_length > max_length:
        pattern = r"(.{%d,})\1{%d,}" % (min_length, min_repeats)

    matches = re.finditer(pattern, sequence)

    tandem_repeats = []

    for m in matches:

        tandem_repeats.append((m.start(), m.end()))

    return tandem_repeats


def collect_all_tandem_repeats(
    dna_sequences: list[str],
    min_length: int = 1,
    max_length: int = -1,
    min_repeats: int = 3,
) -> pd.DataFrame:
    """In each sequence of a given list of DNA sequences, find all tandem_repeats along with:
     - their lenght
     - repeating pattern
     - number of repetitions
     - start and end indices

    Args:
        dna_sequences (list[str]): List of DNA sequences, should contain only A, T, C, G characters.
        min_length (int, optional): The minimum length of the repeating pattern. Defaults to 1.
        max_length (int, optional): The maximum length of the repeating pattern. Defaults to -1,
                                     which indicates no maximum length.
        min_repeats (int, optional): The minimum number of repeats required to consider it a tandem repeat. Defaults to 3.

    Returns:
        pd.DataFrame: A Pandas DataFrame containing the found tandem_repeats and their information, for each sequence.
    """

    all_tandem_repeats = []

    for id, sequence in enumerate(dna_sequences):

        tandem_repeats_idx = find_tandem_repeats_idx(
            sequence, min_length, max_length, min_repeats
        )

        for start, end in tandem_repeats_idx:

            tandem_repeat = sequence[start:end]

            repeating_pattern = find_repeating_pattern(tandem_repeat, min_length)
            n_repetitions = (
                len(tandem_repeat) // len(repeating_pattern) if repeating_pattern else 0
            )

            n_nucleotides = len(repeating_pattern) if repeating_pattern else 0

            all_tandem_repeats.append(
                {
                    "sequence_id": id,
                    "tandem_repeat": tandem_repeat,
                    "length": end - start,
                    "repeating_pattern": repeating_pattern,
                    "n_nucleotides": n_nucleotides,
                    "n_repetitions": n_repetitions,
                    "start": start + 1,  # Adjusting for 1-based index
                    "end": end,
                }
            )

    df_tandem_repeats = pd.DataFrame(all_tandem_repeats)
    df_tandem_repeats.rename_axis("sequence_id", inplace=True)

    return df_tandem_repeats


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

    print("\n~ b. Identify the top 5 most common k-mers (substrings) for k=3, 4, and 5")
    n = 5
    for k in K_TOP_MERS:

        print(f"Compute top {n} {k}-mer count")

        top_kmers = compute_top_kmers_count(dna_sequences, k, n)
        plot_kmer_heatmap(
            top_kmers,
            f"Top {n} {k}-mer Count Heatmap",
            os.path.join(PLOTS_PATH, f"top_{k}_mer_heatmap.jpeg"),
        )

    print("\n~ c. Detect any unusual patterns")
    print("Find palindromes")
    df_palindromes = collect_all_palindromes(dna_sequences, PALINDROMES_MIN_LENGTH)
    print(df_palindromes.head())

    print("\nFind tandem repeats")
    df_tandem_repeats = collect_all_tandem_repeats(
        dna_sequences,
        TANDEM_REPEATS_MIN_LENGTH,
        TANDEM_REPEATS_MAX_LENGTH,
        TANDEM_REPEATS_MIN_REPEATS,
    )

    print("\nWhat are the top 10 longest tandem repeats in the dataset?\n")
    print(df_tandem_repeats.sort_values("length", ascending=False).head(10))

    print(
        "\nWhat are the top 10 repeating patterns most repeated in a single tandem repeat?\n"
    )
    print(df_tandem_repeats.sort_values("n_repetitions", ascending=False).head(10))

    print("\nWhat are the top 10 most common repeating patterns in the dataset?\n")
    print(
        df_tandem_repeats.groupby("repeating_pattern")
        .sum()["n_repetitions"]
        .nlargest(10)
    )
