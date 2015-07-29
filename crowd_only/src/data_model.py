# Tong Shu Li
# First written: 2015-07-02
# Last updated: 2015-07-29
"""
Data models for BioCreative V task 3.

In order to make processing the BioCreative
data easier, these data models were created
in order to simplify the tasks of:
1. Parsing the gold standard
2. Concept highlighting
3. Concept co-occurrence detection
4. Relationship verification
5. Work unit generation
"""
from collections import defaultdict

from file_util import read_file
from lingpipe.split_sentences import split_abstract

def is_MeSH_id(uid):
    return len(uid) == 7 and uid[0] in ["C", "D"]

class Annotation:
    """
    A single mention of a concept in a piece of text.
    Annotation positions are indexed to the abstract.
    """
    def __init__(self, uid, stype, text, start, stop):
        if ":" in uid:
            self.uid_type, self.uid = uid.split(':')
        else:
            self.uid_type = "unknown"
            self.uid = uid

        self.stype = stype.lower()
        assert self.stype in ["chemical", "disease"]

        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop, "Annotation indicies reversed!"
        assert len(text) == self.stop - self.start, "Annotation length mismatch!"

    def __repr__(self):
        """
        Allows for easy printing.
        """
        return "{0}: '{1}'({2}) {3}-{4}".format(self.__class__.__name__,
            self.text, self.stype, self.start, self.stop)

    def __cmp__(self, other):
        """
        When sorting Annotation objects, make sure they
        are in ascending order of starting position.
        """
        if hasattr(other, "start"):
            return self.start.__cmp__(other.start)

class Sentence:
    """
    A single sentence from a paper.
    """
    def __init__(self, pmid, idx, text, start, stop, annotations):
        self.pmid = int(pmid)
        self.uid = "{0}_{1}".format(pmid, idx)
        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop
        assert len(text) == self.stop - self.start

        self.annotations = annotations

    def __repr__(self):
        return "<{0}>: PMID {1} '{2}'({3}-{4})\nAnnotations: {5}\n".format(
            self.__class__.__name__, self.pmid, self.text,
            self.start, self.stop, self.annotations)

class Relation:
    """
    A single chemical-induced disease relationship.
    """
    def __init__(self, drug_id, disease_id):
        assert drug_id != "-1" and disease_id != "-1"
        assert is_MeSH_id(drug_id), "Relation chemical has an improper id!"

        # disease ids can be complex ones joined together
        # represent the disease ids as a set
        disease_id = set(disease_id.split('|'))

        for uid in disease_id:
            assert is_MeSH_id(uid)

        self.drug_id = drug_id
        self.disease_id = disease_id

    def __eq__(self, other):
        """
        Equal if the drug ids match exactly and
        at least one of the disease ids are shared
        between the two relations.

        This is because the gold relations only use
        a pair of single MeSH ids, despite the fact
        that the annotations use complexed MeSH ids.

        WARNING:
            Defining the equals function in this manner
            breaks transitivity. That is, if we have
            three objects A, B, and C, then if
                A == B and B == C, then
                A == C IS NOT TRUE!!

            This is annoying, but the BioCreative data
            is structured badly, so there's nothing
            I can do...
        """
        if isinstance(other, self.__class__):
            return (self.drug_id == other.drug_id
                and len(self.disease_id & other.disease_id) > 0)

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "<{0}>: {1}->{2}".format(
            self.__class__.__name__, self.drug_id, self.disease_id)

class Paper:
    """
    A single academic publication.
    Contains:
        1. The title as a string.
        2. The abstract as a string.
        3. A list of Sentences containing both the title
            and body of the abstract. The first sentence
            is the title.
        4. A list of all chemical and disease annotations
            in the title and abstract sorted in order of
            starting index.
        5. (Optional) A list of gold standard relations.
        6. A set of unique chemical identifiers.
        7. A set of unique disease identifiers.
    """
    def __init__(self, pmid, title, abstract, annotations, relations = []):
        self.pmid = int(pmid)
        self.title = title
        self.abstract = abstract

        self.annotations = sorted(annotations)
        self.relations = relations # may be empty when not parsing gold

        self.chemicals, self.diseases = self.get_unique_concepts(annotations)

        self.sentences = self.split_sentences()

        assert self.has_correct_annotations()

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

    def split_sentences(self):
        """
        Splits the abstract up into individual sentences,
        and determines which concept annotations reside
        within each sentence.

        Time complexity:
            O(N + M) where N is the number of sentences,
            and M is the number of annotations.
        """
        all_sentences = [self.title] + split_abstract(self.abstract)

        full_text = "{0} {1}".format(self.title, self.abstract)

        sent_idx = 0 # starting index of current sentence
        annot_idx = 0 # index of annotation that is within current sentence

        res = []
        for i, sentence in enumerate(all_sentences):
            # The sentence splitter isn't perfect. It recognizes "i.v." as a
            # sentence. Since there can be multiple instances of "sentences"
            # like "i.v." (e.g., PMID 10840460), we need to make sure that
            # we are checking for the first instance starting at the current
            # position (since find always finds the first instance otherwise).
            assert full_text.find(sentence, sent_idx) == sent_idx, "pmid {0} '{1}'".format(
                self.pmid, sentence)

            sent_stop = sent_idx + len(sentence)

            start_annot = annot_idx
            M = len(self.annotations)
            while annot_idx < M and self.annotations[annot_idx].stop <= sent_stop:
                annot_idx += 1

            # should be one past
            res.append(Sentence(self.pmid, i, sentence,
                sent_idx, sent_stop, self.annotations[start_annot : annot_idx]))

            sent_idx += len(sentence) + 1 # all sentences separated by one space

        return res

    def get_work_units(self):
        """
        Returns all possible unique drug-disease combinations
        for this paper.
        """
        return [(drug_id, disease_id) for drug_id in self.chemicals for disease_id in self.diseases]

    def has_correct_annotations(self):
        """
        Checks that the paper's annotations match the
        stated positions in the text.
        """
        text = "{0} {1}".format(self.title, self.abstract)
        for annotation in self.annotations:
            assert text[annotation.start : annotation.stop] == annotation.text

        return True

    def has_relation(self, poss_relation):
        """
        Checks if the provided possible Relationship object
        matches any of the gold standard relationships for this
        paper.

        Note:
            It is not possible to use a set to do the checking
            operation here, because it is not possible to make
            the hashes of two objects the same when they are
            defined by be equal by the overridden equals operator
            for Relation objects.

            This solution is slow, but at least it's correct.
        """
        return poss_relation in self.relations

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
