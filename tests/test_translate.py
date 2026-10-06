"""Tests for skm_tools.translate"""

import gzip

from skm_tools.translate import load_translation_file

TABLE = "ath_source\tstu\nAT1G01010\tSoltu.DM.01G000010\n"


def test_load_translation_file(tmp_path):
    path = tmp_path / "stu.tsv"
    path.write_text(TABLE)
    assert load_translation_file("stu", path).to_dict("records") == [
        {"ath_source": "AT1G01010", "stu": "Soltu.DM.01G000010"}]


def test_load_translation_file_gzipped(tmp_path):
    # as downloaded from skm.nib.si
    path = tmp_path / "stu.tsv.gz"
    path.write_bytes(gzip.compress(TABLE.encode()))
    assert load_translation_file("stu", path)["stu"].tolist() == ["Soltu.DM.01G000010"]
