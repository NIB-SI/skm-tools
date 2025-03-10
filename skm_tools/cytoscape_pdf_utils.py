import io, re

import py4cytoscape as p4c

from pypdf import PdfWriter, PdfReader
from pdfCropMargins import crop
from reportlab.pdfgen import canvas

from .cytoscape_utils import export_network


def _create_text_pdf(text, mb, font, font_size):

    packet = io.BytesIO()

    can = canvas.Canvas(packet, pagesize=(mb.width, mb.height))
    try:
        can.setFont(font, font_size)
    except KeyError as e:
        print(f"Font {font} not available. Using {can._fontname}. ")
        font = can._fontname
        can.setFont(font, font_size)

    text_width = canvas.pdfmetrics._fonts[font].stringWidth(text, size=font_size)
    # rtl_fonts['Helvetica'].stringWidth(text, size=font_size)

    center = (mb.right - mb.left)/2 + mb.left
    text_left = center - text_width/2
    bottom = mb.bottom + font_size/4


    can.drawString(text_left, bottom, text)
    can.save()

    packet.seek(0)
    new_pdf = PdfReader(packet)

    return new_pdf

def export_collection_to_pdfs(collection_suid, folder):
    '''Collection to a pdf per network'''

    networks = sorted(p4c.collections.get_collection_networks(collection_suid))

    i = 0
    if not folder.exists():
        folder.mkdir()

    for network in networks:
        network_name = p4c.get_network_name(network)
        print(network, network_name)

        cropped_pdf = folder / f"{re.sub('[^a-zA-Z0-9]+', '_', network_name).strip('_')}_{network}.pdf"
        pdf = folder / (cropped_pdf.stem + "_cropped" + cropped_pdf.suffix)

        export_network(network, pdf, format="pdf")
        crop(['--noundosave', '-o', str(cropped_pdf), str(pdf)])

        pdf.unlink()

    print(f'Collection saved to {str(folder)}')

def export_collection_to_single_pdf(collection_suid, filename, font_size=20, caption=True, font='Helvetica'):
    '''Collection to single pdf document'''
    '''requires pdf libraries'''

    networks = sorted(p4c.collections.get_collection_networks(collection_suid))

    # TODO rewrite to use "with tempfile.TemporaryDirectory() as tmpdir:"
    i = 0
    while (filename.parent / f"tmp{i}").exists():
        i += 1
    tmp_folder = filename.parent / f"tmp{i}"
    tmp_folder.mkdir()

    writer = PdfWriter()
    for network in networks:
        network_name = p4c.get_network_name(network)
        print(network, network_name)
        pdf = tmp_folder / f"{re.sub('[^a-zA-Z0-9]+', '_', network_name).strip('_')}_{network}.pdf"
        export_network(network, pdf, format="pdf")

        cropped_pdf = tmp_folder / (pdf.stem + "_cropped" + pdf.suffix)
        crop(['--noundosave', '-o', str(cropped_pdf), str(pdf)])

        if caption:
            cropped_pdf2 = tmp_folder / (pdf.stem + "_2nd_cropped" + pdf.suffix)
            crop(['--noundosave', '-a4', '0', f'-{font_size}', '0', '0', '-o', str(cropped_pdf2), str(cropped_pdf)])

            cropped_pdf.unlink()
            cropped_pdf = cropped_pdf2

        reader = PdfReader(cropped_pdf)
        page = reader.pages[0]

        if caption:
            mb = page.cropbox
            caption_pdf = _create_text_pdf(network_name, mb, font, font_size)
            new_page = caption_pdf.pages[0]
            new_page.cropbox = page.cropbox
            page.merge_page(caption_pdf.pages[0])

        writer.add_page(page)

        pdf.unlink()
        cropped_pdf.unlink()

    with open(filename, "wb") as output_stream:
        writer.write(output_stream)
    writer.close()

    tmp_folder.rmdir()

    print(f'Collection save to {filename}')