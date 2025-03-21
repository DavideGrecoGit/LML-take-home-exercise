import os
import pandas as pd
from main import compute_gc_content, plot_gc_content

DNA_SEQUENCES = [
    "AGCTAGCTAG",
    "GGCC",
    "ATATATAT",
    "",
]

EXPECTED_GC_CONTENT = pd.DataFrame(
    [
        {"gc_content": 50.0},  # (G=3 + C=2) / len=10 * 100
        {"gc_content": 100.0},
        {"gc_content": 0.0},
        {"gc_content": 0.0},
    ]
)
EXPECTED_GC_CONTENT.rename_axis("sequence_id", inplace=True)


def test_compute_gc_content():

    result = compute_gc_content(DNA_SEQUENCES)

    pd.testing.assert_frame_equal(result, EXPECTED_GC_CONTENT)


def test_plot_gc_content():

    plot_path = "./test_gc_content.jpg"

    plot_gc_content(EXPECTED_GC_CONTENT, plot_path)

    assert os.path.exists(plot_path)
    os.remove(plot_path)
