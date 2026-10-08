"""Tests for skm_tools.pdf_utils (needs the pdf extra)."""

import pytest

pytest.importorskip("pdfCropMargins")
pypdf = pytest.importorskip("pypdf")
canvas = pytest.importorskip("reportlab.pdfgen.canvas")

from skm_tools.pdf_utils import combine_pdfs, crop_pdf


def _pdf(path, text):
    """A 600 x 400 page with a small text in the middle (lots of margin to crop)."""
    c = canvas.Canvas(str(path), pagesize=(600, 400))
    c.drawString(280, 200, text)
    c.save()
    return path


def _size(path):
    box = pypdf.PdfReader(path).pages[0].cropbox
    return float(box.width), float(box.height)


def test_crop_pdf_in_place(tmp_path):
    path = _pdf(tmp_path / "a.pdf", "network A")
    assert crop_pdf(path) == path
    width, height = _size(path)
    assert width < 300 and height < 100
    assert [p.name for p in tmp_path.iterdir()] == ["a.pdf"]  # no temporary files left


def test_crop_pdf_margins(tmp_path):
    path = _pdf(tmp_path / "a.pdf", "network A")
    _, height = _size(crop_pdf(path, tmp_path / "plain.pdf"))
    _, with_margin = _size(crop_pdf(path, tmp_path / "margin.pdf", margins=(0, 20, 0, 0)))
    assert with_margin == pytest.approx(height + 20, abs=1)
    assert _size(path) == (600, 400)  # the original is unchanged


def test_combine_pdfs_with_captions(tmp_path):
    files = [_pdf(tmp_path / f"{n}.pdf", f"network {n}") for n in "AB"]
    out = combine_pdfs(files, tmp_path / "all.pdf", captions=["JA - Heat", "SA - Heat"])
    reader = pypdf.PdfReader(out)
    assert len(reader.pages) == 2
    assert "JA - Heat" in reader.pages[0].extract_text()
    assert "network B" in reader.pages[1].extract_text()


def test_combine_pdfs_caption_count(tmp_path):
    with pytest.raises(ValueError, match="captions"):
        combine_pdfs([_pdf(tmp_path / "a.pdf", "A")], tmp_path / "all.pdf", captions=[])
