"""
Tong Shu Li
Last updated: 2015-10-19
"""

import os
import pickle
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

def save_file(location, value = None):
    """Saves an object or reads a saved object.

    If provided with just the location, then
    reads the saved object. If an object is provided,
    then that object is saved at the given location.
    """
    if value is None:
        if not os.path.exists(location):
            return None

        with open(location, "rb") as fin:
            res = pickle.load(fin)

        return res

    # save value to file
    with open(location, "wb") as fout:
        pickle.dump(value, fout)
