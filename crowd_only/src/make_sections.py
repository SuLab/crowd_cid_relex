# Tong Shu Li
# Created on Wednesday 2015-07-22
# Last updated 2016-01-08
"""
This program splits the abstract according to the original sections in order
to improve readability.

Section names from the training and development sets are used.
"""
import os

from .lingpipe.file_util import read_file

def create_sections(sentences):
    def has_heading(sentence):
        return any(
            sentence.find("{}:".format(heading)) != -1
            for heading in all_section_names
        )

    fname = os.path.abspath(os.path.join("..", "data", "all_uniq_section_names.txt"))
    all_section_names = [line for line in read_file(fname)]

    positions = {i for i, sentence in enumerate(sentences) if has_heading(sentence)}
    positions |= {0, len(sentences)}
    positions = sorted(list(positions))

    pieces = [
        "<p>{}</p>".format(
            " ".join(sentences[positions[i] : positions[i+1]])
        )
        for i in range(len(positions) - 1)
    ]

    return "".join(pieces)
