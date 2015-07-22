# Tong Shu Li
# Created on Wednesday 2015-07-22
# Last updated 2015-07-22
"""
This program takes the highlighted abstract and
splits it into sections based on the headings
so that the text is easier to read.

Previously this was done by querying the PubMed
API to determine what the section names were,
but now the section names have been cached instead.
"""
from file_util import read_file

def create_sections(abstract):
    # read the list of section names
    fname = "data/all_uniq_section_names.txt"
    all_section_names = [line for line in read_file(fname)]

    positions = [len(abstract)]
    for heading in all_section_names:
        idx = abstract.find("{0}:".format(heading))
        if idx != -1:
            positions.append(idx)

    if len(positions) == 1:
        return abstract

    positions = sorted(positions)

    pieces = []
    for i in range(len(positions) - 1):
        text = abstract[positions[i] : positions[i+1]]
        pieces.append("<p>{0}</p>".format(text))

    return "".join(pieces)
