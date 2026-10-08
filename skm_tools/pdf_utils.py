'''Crop PDF files and combine them into one, with captions (e.g. exported network views).

Requires the ``pdf`` extra (``pip install skm-tools[pdf]``). To export the networks of a
Cytoscape collection as cropped PDFs, see :func:`skm_tools.cytoscape_utils.export_collection`
(``format="PDF", crop=True``).
'''

import io
import logging
import os
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)


def crop_pdf(path, out=None, margins=None):
    '''Crop the margins of a PDF (with pdfCropMargins).

    Parameters
    ----------
    path : str or pathlib.Path
        PDF to crop.
    out : str or pathlib.Path, optional
        File to write (default: replace `path`).
    margins : tuple of float, optional
        Space to add after cropping (left, bottom, right, top), in points, e.g.
        ``(0, 20, 0, 0)`` for a caption below. Default: none.

    Returns
    -------
    pathlib.Path
        The cropped file.
    '''
    from pdfCropMargins import crop

    path = Path(path)
    out = Path(out) if out is not None else path
    args = ["--noundosave"]
    if margins is not None:
        # pdfCropMargins: -a4 adds (negative: removes) space; add = negative "absolute" crop
        args += ["-a4"] + [str(-m) for m in margins]
    fd, tmp = tempfile.mkstemp(dir=out.parent, suffix=".pdf")
    os.close(fd)
    try:
        crop(args + ["-o", tmp, str(path)])
        os.replace(tmp, out)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    return out


def _caption_page(text, box, font, font_size):
    '''A PDF page of the size of `box`, with `text` centred at the bottom.'''
    from pypdf import PdfReader
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfgen import canvas

    packet = io.BytesIO()
    can = canvas.Canvas(packet, pagesize=(float(box.width), float(box.height)))
    try:
        pdfmetrics.getFont(font)
    except KeyError:
        logger.warning("Font %s not available; using %s.", font, can._fontname)
        font = can._fontname
    can.setFont(font, font_size)

    text_width = pdfmetrics.stringWidth(text, font, font_size)
    center = (float(box.right) - float(box.left)) / 2 + float(box.left)
    can.drawString(center - text_width / 2, float(box.bottom) + font_size / 4, text)
    can.save()

    packet.seek(0)
    return PdfReader(packet).pages[0]


def combine_pdfs(files, filename, captions=None, font_size=20, font="Helvetica"):
    '''Combine PDFs into one, a page each (the first page of each file), with captions.

    Parameters
    ----------
    files : iterable of str or pathlib.Path
        PDFs to combine, in order (e.g. from
        :func:`skm_tools.cytoscape_utils.export_collection`). They are not changed.
    filename : str or pathlib.Path
        PDF file to write (overwritten if it exists).
    captions : list of str, optional
        A caption for each page, written below it (default: none).
    font_size : float
        Caption font size (default 20).
    font : str
        Caption font (default Helvetica; the default font if unavailable).

    Returns
    -------
    pathlib.Path
        `filename`.
    '''
    from pypdf import PdfReader, PdfWriter

    files = [Path(f) for f in files]
    if captions is not None and len(captions) != len(files):
        raise ValueError(f"{len(captions)} captions for {len(files)} files.")

    writer = PdfWriter()
    with tempfile.TemporaryDirectory() as tmp:
        for i, path in enumerate(files):
            if captions is not None:
                # room for the caption below
                path = crop_pdf(path, Path(tmp) / f"{i}.pdf", margins=(0, font_size, 0, 0))
            page = PdfReader(path).pages[0]
            if captions is not None:
                page.merge_page(_caption_page(captions[i], page.cropbox, font, font_size))
            writer.add_page(page)

        filename = Path(filename)
        with open(filename, "wb") as handle:
            writer.write(handle)
    writer.close()
    return filename
