"""
Tong Shu Li
First written 2015-07-02
Last updated 2015-08-20

Given a string representing an abstract,
splits it into sentences using LingPipe.
"""
import os
import subprocess

from .file_util import read_file

def split_abstract(abstract):
    """
    Uses two temporary files to talk to the Java
    LingPipe program for sentence splitting.

    Returns a list with the individual sentences.
    """
    orig_dir = os.getcwd()
    real_dir = os.path.dirname(os.path.realpath(__file__))

    os.chdir(real_dir)

    in_fname = "temp_in_file.txt"
    with open(in_fname, "w") as in_file:
        in_file.write(abstract)

    out_fname = "temp_out_file.txt"
    res = subprocess.call(["java", "SplitAbstract", in_fname, out_fname])
    assert res == 0, "Java SplitAbstract did not exit with code 0."

    sentences = [line for line in read_file(out_fname, real_dir)]

    os.chdir(orig_dir)

    return sentences
