# last updated 2015-06-23 tong shu li
"""
This file contains the classes for how we store the information
for the BioCreative V task.
"""

import sys
sys.path.append("/home/toby/Code/util/")
from file_util import read_file

class Annotation:
    """
    A single mention of a concept in a piece of text.
    Annotations are position bound to their abstract.
    """
    def __init__(self, uid, stype, text, start, stop):
        # for use with Dnorm and tmChem:
        if ":" in uid:
            uid = uid.split(':')[1] # get rid of "MESH:", "OMIM:"

        self.uid = uid
        self.stype = stype.lower()
        assert self.stype in ["chemical", "disease"]
        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            # check if positions are the same
            return self.start == other.start and self.stop == other.stop

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def output(self):
        print self.uid
        print self.start
        print self.stop
        print self.text
        print

class Relation:
    """
    A single chemical-induced-disease relationship.
    These relations are bound to their respective
    abstracts.
    """
    def __init__(self, drug, disease):
        assert drug != "-1"
        assert disease != "-1"
        self.drug = drug
        self.disease = disease

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def output(self):
        print self.drug, self.disease

def make_annotations(annotations):
    """
    Annotations with an identifier of -1 are ignored
    because they never show up in a relationship.
    """
    res = []
    for group in annotations:
        if group[5] != "-1":
            res.append(Annotation(group[5], group[4], group[3], group[1], group[2]))

    return res

def make_relations(relations):
    """
    Generate the list of relationships for any
    single abstract.
    """
    return [Relation(group[2], group[3]) for group in relations]

class Paper:
    """
    A single paper containing:
        1. The title and abstract body.
        2. The annotations of concepts.
        3. The relationships between concepts.
    """
    def __init__(self, pmid, title, abstract, annotations, relations):
        self.pmid = pmid
        self.title = title
        self.abstract = abstract
        self.annotations = make_annotations(annotations)
        self.relations = make_relations(relations)

    def output(self):
        print self.pmid
        print len(self.annotations), len(self.relations)

def parse_input(loc, fname):
    """
    Parses the given input file and returns a list
    of Paper objects.
    """
    papers = []

    counter = 0
    annotations = []
    relations = []
    for i, line in enumerate(read_file(fname, loc)):
        if len(line) == 0:
            # time to finish up this paper and prepare a new one
            papers.append(Paper(pmid, title, abstract, annotations, relations))

            counter = 0

            annotations = []
            relations = []
        else:
            if 0 <= counter <= 1:
                vals = line.split('|')
                assert len(vals) == 3
            else:
                vals = line.split('\t')

            if counter == 0:
                assert vals[1] == 't'
                pmid = vals[0]
                title = vals[2]
            elif counter == 1:
                assert vals[1] == 'a'
                abstract = vals[2]
            elif len(vals) == 4:
                relations.append(vals)
            else:
                assert 6 <= len(vals) <= 7, i
                annotations.append(vals) # 6 or 7 fields

            counter += 1

    return papers
