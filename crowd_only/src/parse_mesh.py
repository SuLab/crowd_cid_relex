"""
Tong Shu Li
Created on: 2015-09-01
Last updated: 2015-10-14

Parses the MeSH XML file into a useable format.
"""
from collections import defaultdict
import os
import pickle
import xml.etree.cElementTree as ET

def parse_mesh_xml(fname):
    # all tree numbers for each concept are unique
    tree = ET.parse(fname)
    root = tree.getroot()

    hierarchy = defaultdict(list)
    concept_name = dict()
    for record in root:
        mesh_id = record.find("DescriptorUI").text

        name_root = record.find("DescriptorName")
        assert len(name_root) == 1
        name = name_root.find("String").text

        if mesh_id not in concept_name:
            concept_name[mesh_id] = name
        else:
            assert concept_name[mesh_id] == name

        tree_num_root = record.find("TreeNumberList")
        if tree_num_root is None:
            hierarchy[mesh_id] = set()
        else:
            for tree_num in record.find("TreeNumberList"):
                hierarchy[mesh_id].append(tree_num.text)

    return (concept_name, hierarchy)

def parse_supplement(fname):
    """Determine the official names of the supplemental concepts."""
    tree = ET.parse(fname)
    root = tree.getroot()

    mapping = dict()
    for record in root.iter("SupplementalRecord"):
        uid = record.find("SupplementalRecordUI").text
        name = record.find("./SupplementalRecordName/String").text

        mapping[uid] = name

    return mapping

def load_mesh(fname):
    assert fname in ["supp", "hierarchy"], "Invalid MeSH pickle name"

    cur = os.path.dirname(__file__)

    loc = os.path.join(cur, "..", "data", "mesh_ontology", "mesh_{}.pickle".format(fname))
    with open(os.path.abspath(loc), "rb") as fin:
        val = pickle.load(fin)

    return val

def main():
    loc = os.path.join("..", "data", "mesh_ontology")

    fname = os.path.abspath(os.path.join(loc, "supp2015.xml"))
    supplement = parse_supplement(fname)

    with open(os.path.join(loc, "mesh_supp.pickle"), "wb") as fout:
        pickle.dump(supplement, fout)

#-------------------------------------------------------------------------------

    fname = os.path.abspath(os.path.join(loc, "mesh2015.xml"))
    concept_name, hierarchy = parse_mesh_xml(fname)

    with open(os.path.join(loc, "mesh_hierarchy.pickle"), "wb") as fout:
        pickle.dump((concept_name, hierarchy), fout)

if __name__ == "__main__":
    main()
