"""
Tong Shu Li
2015-07-02

Given a string representing an abstract,
splits it into sentences using LingPipe.
"""
import subprocess

import sys
sys.path.append("/home/toby/Code/util/")
from file_util import read_file

def split_abstract(abstract):
    """
    Uses two temporary files to talk to the Java
    LingPipe program for sentence splitting.

    Returns a list with the individual sentences.
    """
    in_fname = "temp_in_file.txt"
    with open(in_fname, "w") as in_file:
        in_file.write(abstract)

    out_fname = "temp_out_file.txt"
    res = subprocess.call(["java", "SplitAbstract", in_fname, out_fname])
    assert res == 0, "Java SplitAbstract did not exit with code 0."

    return [line for line in read_file(out_fname)]
