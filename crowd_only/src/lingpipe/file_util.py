"""
Tong Shu Li
Last updated: 2015-08-20
"""

import os
import shutil

def read_file(file_name, file_loc = os.getcwd()):
    with open(os.path.join(file_loc, file_name), "r") as file:
        for line in file:
            yield line.rstrip('\n')

def make_dir(loc):
    if not os.path.exists(loc):
        os.makedirs(loc)

def move_and_replace(file, dest):
    if os.path.exists(dest + file):
        os.remove(dest + file)
    shutil.move(file, dest)

def exists(file_name, file_loc = os.getcwd()):
    return os.path.exists(os.path.join(file_loc, file_name))
