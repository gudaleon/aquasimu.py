"""The listing reader is exercised on a file this test writes. No AQUASIM output is stored."""
from aquasimu.io.listing import read_listing


def test_read_listing_parses_column_pairs(tmp_path):
    path = tmp_path / "listing.txt"
    path.write_text(
        "Legend:\t\tA\t\tB\n"
        "\t0\t1\t0\t2\n"
        "\t1\t3\t1\t4\n"
        "\n"
    )
    legend, curves = read_listing(path)
    assert list(legend) == ["A", "B"]
    assert curves[0][0].tolist() == [0.0, 1.0]
    assert curves[0][1].tolist() == [1.0, 3.0]
    assert curves[1][1].tolist() == [2.0, 4.0]
