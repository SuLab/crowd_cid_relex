"""
Tong Shu Li
Created on: 2015-10-05
Last updated: 2015-10-15

Use the MeSH hierarchy to determine which concepts are more specific.
"""
from collections import namedtuple
from itertools import groupby
import os

from .data_model import Ontology_ID
from .data_model import Simple_Rel
from .parse_mesh import load_mesh

concept_name, hierarchy = load_mesh("hierarchy")

def mesh_filter(concepts):
    """Removes redundant disease concepts from a list using the MeSH hierarchy.

    Given a list of MeSH disease IDs, removes those which are strict parents
    (prefixes) of other concepts in the list. Returns a set of ids.

    N^2 runtime on the number of concepts to verify.
    """
    def redundant(tree_num):
        for label in seen:
            if label.startswith(tree_num):
                return True

        return False

    Ref = namedtuple("Ref", ["tree_num", "mesh_id"])

    tree_nums = []
    for concept in concepts:
        tree_nums += [Ref(tree_num, concept) for tree_num in hierarchy[concept]]

    tree_nums = sorted(tree_nums, key = lambda val: len(val.tree_num), reverse = True)

    # grab each identifier one by one
    seen = set()
    ans = set()
    for tree_num, concept in tree_nums:
        if concept in ans:
            seen.add(tree_num)
        elif not redundant(tree_num):
            ans.add(concept)
            seen.add(tree_num)

    return ans

def filter_relations(relations):
    """Filter redundant relationships from a list of Simple_Rel objects.

    Using the MeSH ontology, redundant relations between one chemical and
    multiple diseases are removed. A relation is considered redundant if the
    disease is a more general than another disease related to the same chemical.

    Example:
        Kidney failure, chronic (D007676) is more specific than Kidney diseases
        (D007674). If we have the relations:

        [(A, D007676), (A, D007674), (B, D007674)]

        the output will be:

        [(A, D007676), (B, D007674)]

        Since D007674 is more general, but no more specific relation to B
        exists, and is therefore only removed when related to A.
    """

    relations = sorted(relations, key = lambda val: val.chemical.flat_repr)

    # assume all the same pmid
    pmid = relations[0].pmid

    ans = []
    for chemical, group in groupby(relations, lambda val: val.chemical):
        diseases = mesh_filter([rel.disease.uid for rel in group])
        ans += [Simple_Rel(pmid, chemical, Ontology_ID(disease)) for disease in diseases]

    return ans
