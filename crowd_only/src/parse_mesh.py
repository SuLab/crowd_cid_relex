"""
Tong Shu Li
Created on: 2015-09-01
Last updated: 2015-10-05

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

def main():
    parent = os.path.abspath("..")
    fname = os.path.join(parent, "data/mesh_ontology/mesh2015.xml")

    concept_name, hierarchy = parse_mesh_xml(fname)

    # save to file
    data_dir = os.path.join(parent, "data/mesh_ontology")
    with open(os.path.join(data_dir, "mesh_names.pickle"), "wb") as fout:
        pickle.dump(concept_name, fout)

    with open(os.path.join(data_dir, "mesh_hierarchy.pickle"), "wb") as fout:
        pickle.dump(hierarchy, fout)

if __name__ == "__main__":
    main()
