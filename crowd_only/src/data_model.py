# Tong Shu Li
# First written: 2015-07-02
# Last updated: 2015-12-16
"""Data models for BioCreative V chemical-disease relation task.

These data models were created to simplify the tasks of:
1. Parsing the gold standard
2. Concept co-occurrence detection
3. Relationship verification
4. Work unit generation
"""
from collections import defaultdict
from itertools import islice

from .lingpipe.file_util import read_file
from .lingpipe.file_util import save_file
from .lingpipe.split_sentences import split_abstract

def is_MeSH_id(uid):
    return len(uid) == 7 and uid[0] in ["C", "D"]

#-------------------------------------------------------------------------------

class Base:
    def __init__(self, uid):
        self.uid = uid

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__

        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not(self == other)

        return NotImplemented


class OntologyID(Base):
    """A single identifier from an existing biomedical ontology."""

    def __init__(self, text):
        """Create an OntologyID from a string representation."""

        assert text.count(":") <= 1, "ID {} is misformatted!".format(text)
        assert "|" not in text

        if ":" in text:
            self.uid_type, self.uid = text.split(":")
            if self.uid_type == "MESH":
                assert is_MeSH_id(self.uid), "{} not MeSH!".format(text)
        else:
            self.uid = text
            self.uid_type = "MESH" if is_MeSH_id(self.uid) else "unknown"

        self.flat_repr = "{}:{}".format(self.uid_type, self.uid)

    def __repr__(self):
        return "<{0}>: {1}".format(self.__class__.__name__, self.flat_repr)

    def __hash__(self):
        return hash(self.flat_repr)

    def __lt__(self, other):
        """Sort OntologyIDs by their string representation."""
        if isinstance(other, self.__class__):
            return self.flat_repr < other.flat_repr

        return NotImplemented


class MultiID(Base):
    """One or more OntologyIDs used to identify an Annotation.

    The BioCreative dataset annoyingly can assign multiple MeSH or other
    identifiers to a single annotation. However, in the gold standard relations
    only single identifiers are used for concepts. Therefore some sort of
    expansion is needed when converting annotations to relations which can be
    compared with the gold standard.

    This class aims to store all the IDs for an annotation and give easy access
    to the MeSH ids contained within.

    E.g., PMID 8638876 (development set)
        The annotation used for a relation is "D002544|-1". No other
        annotations contain D002544 as an identifier.
    """

    def __init__(self, text):
        """Create a Multi_ID from a string."""
        self.uid = tuple(sorted(OntologyID(v) for v in text.split("|")))
        self.flat_repr = "|".join(v.flat_repr for v in self.uid)

    def get_mesh_only(self):
        """Grabs the MeSH OntologyIDs only."""
        return tuple(v for v in self.uid if v.uid_type == "MESH")

    def has_mesh(self):
        """Checks if any of the OntologyIDs are MeSH."""
        return any(v.uid_type == "MESH" for v in self.uid)

    def update_ids(self, new_id):
        """Used to update an annotation's identity during acronym resolution."""
        assert isinstance(new_id, MultiID)
        self.uid = new_id.uid
        self.flat_repr = new_id.flat_repr

    def __hash__(self):
        return hash(self.uid)

    def __repr__(self):
        return "<{}>: {}".format(self.__class__.__name__, self.flat_repr)


class Position(Base):
    def __init__(self, text, start, stop):
        self.text = text
        self.start = int(start)
        self.stop = int(stop)

        assert self.start < self.stop, "{0} indicies reversed!".format(self)
        assert len(text) == self.stop - self.start, "{0} length mismatch!".format(self)

    def __lt__(self, other):
        """Sort by position."""
        if isinstance(other, self.__class__):
            return self.start < other.start or (self.start == other.start
                and self.stop < other.stop)

        return NotImplemented


class Annotation(Position):
    """A single mention of a concept in a piece of text."""

    def __init__(self, uid, stype, text, start, stop):
        Base.__init__(self, MultiID(uid))
        Position.__init__(self, text, start, stop)

        self.stype = stype.lower()
        assert self.stype in ["chemical", "disease"], "Bad semtype: {}".format(self)

    def __repr__(self):
        return "<{}>: '{}'({}-{}, {})".format(self.__class__.__name__,
            self.text, self.start, self.stop, self.uid)

    def __hash__(self):
        return hash((self.text, self.start, self.stop, self.stype, self.uid))


class Sentence(Position):
    """A single sentence.

    All unique chemical-disease relationships between annotations contained by
    this sentence are classified into CID and non-CID relations.

    This sentence's non-CID relations minus the CID relations true at the
    abstract level need to be verified in a sentence-level task.
    """
    def __init__(self, pmid, idx, text, start, stop, annotations):
        Base.__init__(self, "{}_{}".format(pmid, idx))
        self.pmid = int(pmid)

        # start = position of sentence in the abstract
        Position.__init__(self, text, start, stop)

        # a list of the concept annotations within this sentence
        # annotations should already be sorted at the paper stage
        self.annotations = annotations

        # generate the list of CID and non-CID relations bound to this sentence
        self.poss_relations = self.classify_relations()

    def __repr__(self):
        return "<{}>: PMID:{} '{}'({}-{})\nAnnotations: {}\n".format(
            self.__class__.__name__, self.pmid, self.text,
            self.start, self.stop, self.annotations)

    def is_CID_relation(self, chemical, disease):
        """Determine if a chemical annotation and a disease annotation follow
        the CID structure.
        """
        return ((chemical.stop < disease.start)
            and (disease.start - chemical.stop <= 15)
            and ("induce" in self.text[chemical.stop - self.start :
                disease.start - self.start].lower())
        )

    def classify_relations(self):
        """Classify all unique chemical-disease identifier pairs in this
        sentence as CID or non-CID relations.

        The CID and non-CID relation identifier pairs are mutually exclusive.
        Only MeSH identifiers will be generated, since the gold standard only
        lists relations between MeSH concepts. Any non-MeSH concepts can be
        skipped.
        """
        all_relations = defaultdict(set)
        for annot_A in self.annotations:
            if annot_A.stype == "chemical" and annot_A.uid.has_mesh():
                for annot_B in self.annotations:
                    if annot_B.stype == "disease" and annot_B.uid.has_mesh():
                        rel_type = self.is_CID_relation(annot_A, annot_B)
                        all_relations[rel_type].add((annot_A.uid, annot_B.uid))

        """
        In cases where we have a sentence with the following annotations:
        D C D (where D = disease and C = chemical), we see that the first
        instance of C and D is not in the CID format, and will get added
        to the non-CID set. However, the second instance is in the CID
        format, and gets added to the CID set. To make sure the two sets
        are mutually exclusive, we need to subtract the CID set from the
        larger non-CID set.
        """
        all_relations[False] -= all_relations[True]
        assert all_relations[False].isdisjoint(all_relations[True])
        return all_relations


class SingleRelation(Base):
    """A single CID relation between two single MeSH OntologyIDs.

    Relations contain information about the likely annotation they derived from,
    as well as whether the concepts cooccurred within a sentence.

    Mainly used to represent one gold standard relation. Relations which are
    generated in a forward process from annotations (which may have Multi_IDs)
    should not be kept in a Relation object, but should rather be processed
    into a set of Relation objects.
    """
    def __init__(self, pmid, chem_id, dise_id, sent_cooccur, chem_annot, dise_annot):
        """Create a new Relation object.

        chem_id and dise_id are the single OntologyIDs of the relation (gold
        standard format).

        sent_cooccur is a boolean representing whether the two concepts ever
        cooccurred within any sentences.

        chem_annot = the Multi_ID representing the annotation the chem_id derives
        from

        dise_annot = the Multi_ID representing the annotation the dise_id derives
        from
        """
        self.pmid = int(pmid)

        assert isinstance(chem_id, OntologyID)
        assert isinstance(dise_id, OntologyID)

        assert chem_id.uid_type == "MESH" and dise_id.uid_type == "MESH"

        assert isinstance(chem_annot, MultiID)
        assert isinstance(dise_annot, MultiID)

        assert isinstance(sent_cooccur, bool)

        self.chem = chem_id
        self.dise = dise_id

        self.sent_cooccur = sent_cooccur
        self.chem_orig = chem_annot
        self.dise_orig = dise_annot

    def __repr__(self):
        return "<{}>: {}&{}(PMID:{}, {}, {}, {})".format(self.__class__.__name__,
            self.chem, self.dise, self.pmid,
            int(self.sent_cooccur), self.chem_orig, self.dise_orig)

    def __hash__(self):
        return hash((self.pmid, self.chem, self.dise, self.sent_cooccur,
            self.chem_orig, self.dise_orig))


class RelationGroup:
    """A collection of SingleRelation objects.

    A data structure for holding the set of all gold standard relations
    associated with a paper is needed. In addition, there needs to be a way to
    store all the SingleRelations associated with the crowd predictions.

    This data structure aims to fulfill this need, and is responsible for
    determining the original annotations used for each gold standard relation.
    """
    def __init__(self, lel):
        self.stuff = lel


class Paper:
    """A single academic abstract.

    Contains:
        1. The PubMed identifier as an integer.
        2. The title as a string.
        3. The abstract as a string.
        4. A list of all chemical and disease annotations in the title and
            abstract sorted in increasing order of starting index.
        5. A potentially empty list of gold standard CID relations.
        6. A set of unique chemical identifiers.
        7. A set of unique disease identifiers.
        8. A list of Sentences containing both the title and body of the
            abstract. The 0th sentence is the title. Each sentence contains the
            annotations and relations constrained to that particular sentence.
        9. A set of all the potential chemical-disease relations grouped into
            three mutually exclusive categories:
            - CID relations
            - Non-CID sentence-bound relations
            - Non-sentence bound relations (by definition not CID relations)

            The sum of relations in all three groups should equal the number of
            unique chemical IDs times the number of unique disease IDs.
        10. A function for resolving acronyms.
    """
    def __init__(self, pmid, title, abstract, annotations,
        gold_relations = [], fix_acronyms = False):

        self.pmid = int(pmid)
        self.title = title
        self.abstract = abstract

        self.annotations = sorted(annotations)
        assert self.has_correct_annotations()
        if fix_acronyms:
            self.resolve_acronyms()

        self.chemicals, self.diseases = self.get_unique_concepts()

        # split sentences and generate sentence-bound relations
        self.sentences = self.split_sentences()
        self.poss_relations = self.classify_relations()

        self.gold_relations = gold_relations # may be empty when not parsing gold

    def __repr__(self):
        return ("<{0}>: PMID {1}. {2} annotations, {3} gold relations\n"
            "{4} unique chemical ids, {5} unique disease ids\n"
            "{6} sentences".format(self.__class__.__name__,
            self.pmid, len(self.annotations), len(self.gold_relations),
            len(self.chemicals), len(self.diseases), len(self.sentences))
        )

    def has_correct_annotations(self):
        """Check that the annotation indicies produce the correct text snippet.

        Also check that annotations do not overlap with one another.
        """
        text = "{} {}".format(self.title, self.abstract)
        for annot in self.annotations:
            assert text[annot.start : annot.stop] == annot.text, (
                "Annotation {} text mismatch".format(annot))

        for i, annot in enumerate(self.annotations[:-1]):
            assert self.annotations[i+1].start > annot.stop, (
                "Annotation overlap {}".format(annot))

        return True

    def resolve_acronyms(self):
        """Identifies annotations likely to be acronyms for concepts stated
        earlier in the text and assigns them the same identifier as the original
        concept definition.

        tmChem identifies many text spans as likely chemicals, but is sometimes
        unable to assign a MeSH identifier to the annotation. These annotations
        are often acronyms defined earlier in the text. Assigning the acronyms
        the same MeSH identifier as the original definition improves tmChem
        performance.
        """
        used = [False] * len(self.annotations)
        full_text = "{} {}".format(self.title, self.abstract)

        # if an abbreviation is included in parentheses, then it should
        # follow the definition annotation immediately
        for i, definition in enumerate(self.annotations[ : -1]):
            if not used[i] and definition.uid.has_mesh():
                acronym = self.annotations[i + 1]

                if (acronym.stype == definition.stype
                    and acronym.start == definition.stop + 2
                    and full_text[acronym.start - 1] == "("
                    and full_text[acronym.stop] == ")"):

                    # found an acronym definition

                    used[i] = True
                    for j, annot in enumerate(islice(self.annotations, i + 1, None)):
                        if (annot.stype == definition.stype
                            and not used[i + 1 + j]
                            and not annot.uid.has_mesh()
                            and annot.text == acronym.text):

                            self.annotations[i + 1 + j].update_uid(definition.uid)
                            used[i + 1 + j] = True

    def get_unique_concepts(self):
        """Determine the unique chemical and disease identifiers for this paper.

        Ignores any annotations which are not assigned MeSH identifiers.
        """
        res = defaultdict(set)
        for annotation in self.annotations:
            if annotation.uid.has_mesh():
                res[annotation.stype].add(annotation.uid)

        return (res["chemical"], res["disease"])

    def get_all_possible_relations(self):
        """Return all possible unique chemical-disease identifier pairs."""
        return {
            (chemical_id, disease_id)
            for chemical_id in self.chemicals
                for disease_id in self.diseases
        }

    def split_sentences(self):
        """Split the abstract into individual sentences, and determine which
        concept annotations reside within each sentence.

        Time complexity: O(N + M) where N is the number of sentences and M is
        the number of annotations.
        """
        all_sentences = [self.title] + split_abstract(self.abstract)

        full_text = "{} {}".format(self.title, self.abstract)

        sent_idx = 0 # starting index of current sentence
        annot_idx = 0 # index of annotation that is within current sentence

        res = []
        M = len(self.annotations)
        for i, sentence in enumerate(all_sentences):
            # The sentence splitter isn't perfect. It recognizes "i.v." as a
            # sentence. Since there can be multiple instances of "sentences"
            # like "i.v." (e.g., PMID 10840460), we need to make sure that
            # we are checking for the first instance starting at the current
            # position (since find always finds the first instance otherwise).
            assert full_text.find(sentence, sent_idx) == sent_idx, (
                "PMID {0} {1} text mismatch!".format(self.pmid, sentence))

            sent_stop = sent_idx + len(sentence)

            start_annot = annot_idx
            while annot_idx < M and self.annotations[annot_idx].stop <= sent_stop:
                annot_idx += 1

            # should be one past
            res.append(Sentence(self.pmid, i, sentence,
                sent_idx, sent_stop, self.annotations[start_annot : annot_idx]))

            sent_idx += len(sentence) + 1 # all sentences separated by one space

        return res

    def classify_relations(self):
        """Classify all unique chemical-disease identifier pairs in this
        abstract into three groups:

        1. CID relations, which are sentence bound
        2. Non-CID, sentence-bound relations
        3. Relations which are not sentence bound
        """
        all_rels = self.get_all_possible_relations()

        cid_rels = set()
        sentence_rels = set()
        for sentence in self.sentences:
            cid_rels |= sentence.poss_relations[True]
            sentence_rels |= sentence.poss_relations[False]

        sentence_rels -= cid_rels
        abs_rels = all_rels - cid_rels - sentence_rels

        assert cid_rels.isdisjoint(sentence_rels)
        assert cid_rels.isdisjoint(abs_rels)
        assert abs_rels.isdisjoint(sentence_rels)

        assert all_rels == cid_rels | sentence_rels | abs_rels

        return {
            "CID": cid_rels,
            "sentence_non_CID": sentence_rels,
            "not_sentence_bound": abs_rels
        }


def parse_input(loc, fname, fix_acronyms = True):
    """Parse a PubTator formatted file and return a dict of Paper objects."""

    papers = dict()
    counter = 0
    annotations = []
    relations = []
    for i, line in enumerate(read_file(fname, loc)):
        if not line:
            papers[pmid] = Paper(pmid, title, abstract, annotations,
                relations, fix_acronyms = fix_acronyms)

            counter = -1
            annotations = []
            relations = []
        elif counter < 2:
            vals = line.split('|')
            assert len(vals) == 3, "Bad format for line {}".format(i+1)
            assert vals[1] == ["t", "a"][counter]

            if counter == 0:
                pmid = int(vals[0])
                title = vals[2]
            else:
                assert pmid == int(vals[0])
                abstract = vals[2]
        else:
            vals = line.split('\t')
            assert pmid == int(vals[0])
            if vals[1] == "CID":
                relations.append((OntologyID(vals[2]), OntologyID(vals[3])))
            else:
                # an annotation
                if len(vals) == 5: # no identifier was assigned
                    vals.append("-1")

                assert 6 <= len(vals) <= 7, "Error on line {0}".format(i+1)
                annotations.append(Annotation(vals[5], vals[4], vals[3], vals[1], vals[2]))

        counter += 1

    return papers


def parse_file(save_loc, **kwargs):
    """Uses a cached version of the save file if possible."""
    res = save_file(save_loc)
    if res is not None:
        return res

    res = parse_input(kwargs["loc"], kwargs["fname"],
        fix_acronyms = kwargs["fix_acronyms"])

    save_file(save_loc, res)
    return res
