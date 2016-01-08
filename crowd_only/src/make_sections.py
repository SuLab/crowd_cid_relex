# Tong Shu Li
# Created on Wednesday 2015-07-22
# Last updated 2016-01-08
"""
This program takes the highlighted abstract and
splits it into sections based on the headings
so that the text is easier to read.

Previously this was done by querying the PubMed
API to determine what the section names were,
but now the section names have been cached instead.
"""
from collections import defaultdict
import os

from .lingpipe.file_util import read_file

def create_sections(abstract):
    """
    Some section headings are suffixes of other headings (e.g., METHODS:,
    DESIGN/METHODS:). When this occurs, both headings will be found, but we want
    to take the longer heading.
    """
    fname = os.path.abspath(os.path.join("..", "data", "all_uniq_section_names.txt"))
    all_section_names = [line for line in read_file(fname)]

    earliest_idx = defaultdict(set)
    for heading in all_section_names:
        # assumes that each section name will only ever appear at most once
        start = abstract.find("{}:".format(heading))
        if start != -1:
            # record all headings ending at this position
            stop = start + len(heading)
            earliest_idx[stop].add(start)

    positions = {min(starts) for end_pos, starts in earliest_idx.items()}
    if not positions:
        return abstract

    positions |= {0, len(abstract)}
    positions = sorted(list(positions))

    pieces = [
        "<p>{}</p>".format(abstract[positions[i] : positions[i+1]])
        for i in range(len(positions) - 1)
    ]

    return "".join(pieces)
