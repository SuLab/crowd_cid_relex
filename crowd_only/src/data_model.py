# Tong Shu Li
# First written: 2015-07-02
# Last updated: 2015-07-02
"""
Revised data models for BioCreative V.
"""
from collections import defaultdict

import sys
sys.path.append("/home/toby/Code/util/")
from file_util import read_file

from lingpipe.split_sentences import split_abstract

def is_MeSH_id(uid):
    return len(uid) == 7 and uid[0] in ["C", "D"]

class Annotation:
    """
    A single mention of a concept in a piece of text.
    Annotation positions are bound to the abstract.
    """
    def __init__(self, uid, stype, text, start, stop):
        if ":" in uid:
            self.uid_type, self.uid = uid.split(':')
        else:
            self.uid = uid
            self.uid_type = "unknown"

        self.stype = stype.lower()
        assert self.stype in ["chemical", "disease"]

        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop
        assert len(text) == self.stop - self.start

    def output():
        print "Id: {0}, {1}, {2}, {3}, {4}".format(self.uid, self.stype, self.text, self.start, self.stop)

class Sentence:
    """
    A single sentence from a paper.
    """
    def __init__(self, pmid, idx, text, start, stop, paper_annotations):
        self.pmid = int(pmid)
        self.uid = "{0}_{1}".format(pmid, idx)
        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop
        assert len(text) == self.stop - self.start

        self.annotations = self.find_annotations(paper_annotations)

    def find_annotations(self, annotations):
        """
        Given a list of the paper's annotations, determines
        which are within this sentence.
        """
        res = []
        for annotation in annotations:
            if annotation.start >= self.start and annotation.stop <= self.stop:
                res.append(annotation)

        return res

class Relation:
    """
    A single chemical-induced disease relationship.
    """
    def __init__(self, drug_id, disease_id):
        assert drug_id != "-1" and disease_id != "-1"
        assert is_MeSH_id(drug_id) and is_MeSH_id(disease_id)
        self.drug_id = drug_id
        self.disease_id = disease_id

class Paper:
    """
    A single academic publication.
    Contains:
        1. The title as a string.
        2. The abstract as a string.
        3. A list of Sentences representing the text of
            the body of the abstract.
        4. A list of annotations.
        5. (Optional) A list of gold standard relations.
        6. A set of unique chemical identifiers.
        7. A set of unique disease identifiers.
    """
    def __init__(self, pmid, title, abstract, annotations, relations):
        self.pmid = int(pmid)
        self.title = title
        self.abstract = abstract
        self.annotations = annotations
        self.relations = relations

        self.chemicals, self.diseases = self.get_unique_concepts(annotations)
        self.sentences = self.split_sentences(abstract)

    def get_unique_concepts(self, annotations):
        """
        Determines the unique identifiers of chemicals and diseases
        belonging to this paper.

        Ignores any annotations with an identifier of -1.
        """
        res = defaultdict(set)
        for annotation in annotations:
            if annotation.uid != "-1":
                res[annotation.stype].add(annotation.uid)

        return (res["chemical"], res["disease"])

    def split_sentences(self, abstract):
        """
        Splits the paper's abstract up into individual
        sentences, and determines which annotations are
        in each sentence.
        """
        sentences = split_abstract(abstract)
        full_text = "{0} {1}".format(self.title, self.abstract)

        res = []
        for i, sentence in enumerate(sentences):
            idx = full_text.index(sentence)
            res.append(Sentence(self.pmid, i, sentence, idx, idx + len(sentence), self.annotations))

        return res

def parse_input(loc, fname, is_gold = True):
    """
    Reads a given file and returns a list of Paper
    objects.
    """
    papers = []
    counter = 0
    annotations = []
    relations = []
    for i, line in enumerate(read_file(fname, loc)):
        if len(line) == 0:
            # time to create the paper object
            papers.append(Paper(pmid, title, abstract, annotations, relations))

            counter = 0
            annotations = []
            relations = []
        else:
            if 0 <= counter <= 1:
                vals = line.split('|')
                assert len(vals) == 3, "Title or abstract on line {i} is messed up!".format(i + 1)
            else:
                vals = line.split('\t')

            if counter == 0:
                assert vals[1] == "t", i+1
                pmid = vals[0]
                title = vals[2]
            elif counter == 1:
                assert vals[1] == "a"
                assert vals[0] == pmid
                abstract = vals[2]
            elif is_gold and len(vals) == 4:
                assert vals[1] == "CID"
                relations.append(Relation(vals[2], vals[3]))
            else:
                assert 6 <= len(vals) <= 7
                annotations.append(Annotation(vals[5], vals[4], vals[3], vals[1], vals[2]))

            counter += 1

    return papers
