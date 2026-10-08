"""Tests for skm_tools.utils"""

import gzip
import io
from email.message import Message
from unittest.mock import patch

import pytest

from skm_tools.utils import as_list, download, merge_values


@pytest.mark.parametrize("x, expected", [
    ("ath", ["ath"]), (("ath", "stu"), ["ath", "stu"]), ({"ath"}, ["ath"]), (range(2), [0, 1]),
    (3, [3]), (None, None), ([], []),
])
def test_as_list(x, expected):
    assert as_list(x) == expected


@pytest.mark.parametrize("attribute, values, expected", [
    ("interaction", ["positive-influence", "positive-influence"], ("positive-influence", False)),
    ("interaction", ["positive-influence", "negative-influence"], ("unknown-influence", True)),
    ("directed", [False, True], (True, True)),
    ("rank", [3, None, 1], (1, True)),
    ("reaction_type", ["binding", "catalysis"], ("binding", True)),
    ("reaction_type", [None, "catalysis"], ("catalysis", False)),
    ("forms", [["a"], ["a"]], (["a"], False)),
    ("forms", [None, None], (None, False)),
])
def test_merge_values(attribute, values, expected):
    assert merge_values(attribute, values) == expected


class _Response(io.BytesIO):
    def __init__(self, data, content_type):
        super().__init__(data)
        self.headers = Message()
        self.headers["Content-Type"] = content_type

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


TEXT = b"a\tb\n1\t2\n"


@pytest.mark.parametrize("content_type, data, name", [
    ("application/gzip", gzip.compress(TEXT), "f.tsv.gz"),
    ("application/gzip", gzip.compress(TEXT), "f.tsv"),
    ("text/tab-separated-values", TEXT, "f.tsv"),
    ("text/tab-separated-values", TEXT, "f.tsv.gz"),
])
def test_download_compresses_as_named(tmp_path, content_type, data, name):
    with patch("skm_tools.utils.urlopen", return_value=_Response(data, content_type)):
        path = download("https://example.org/f", tmp_path / "sub" / name)
    content = path.read_bytes()
    assert (gzip.decompress(content) if name.endswith(".gz") else content) == TEXT
    assert [p.name for p in path.parent.iterdir()] == [name]


def test_interrupted_download_leaves_no_file(tmp_path):
    class Broken(_Response):
        def read(self, *args):
            raise ConnectionError("interrupted")

    with patch("skm_tools.utils.urlopen", return_value=Broken(TEXT, "text/plain")):
        with pytest.raises(ConnectionError):
            download("https://example.org/f", tmp_path / "f.tsv")
    assert list(tmp_path.iterdir()) == []
