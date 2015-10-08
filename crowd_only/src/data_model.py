# Tong Shu Li
# First written: 2015-07-02
# Last updated: 2015-10-05
"""Data models for BioCreative V task 3.

These data models were created to simplify the tasks of:
1. Parsing the gold standard
2. Concept co-occurrence detection
3. Relationship verification
4. Work unit generation
"""
from collections import defaultdict
from itertools import islice

from .lingpipe.file_util import read_file
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
            return not self == other

        return NotImplemented

class Ontology_ID(Base):
    """A single identifier from an existing biomedical ontology."""

    def __init__(self, text):
        """Create an Ontology_ID from a string representation."""

        assert text.count(":") <= 1, "ID {0} is misformatted!".format(text)
        assert "|" not in text, "Was given compound identifier!"
        if ":" in text:
            self.uid_type, self.uid = text.split(":")
            if self.uid_type == "MESH":
                assert is_MeSH_id(self.uid), "MeSH ID mismatch for {0}".format(text)
        else:
            self.uid = text
            self.uid_type = "MESH" if is_MeSH_id(self.uid) else "unknown"

        self.flat_repr = "{0}:{1}".format(self.uid_type, self.uid)

    def __repr__(self):
        return "<{0}>: {1}".format(self.__class__.__name__, self.flat_repr)

    def __hash__(self):
        return hash(self.flat_repr)

class Position(Base):
    def __init__(self, text, start, stop):
        self.text = text
        self.start = int(start)
        self.stop = int(stop)

        assert self.start < self.stop, "{0} indicies reversed!".format(self)
        assert len(text) == self.stop - self.start, "{0} length mismatch!".format(self)

class Annotation(Position):
    """A single mention of a concept in a piece of text."""

    def __init__(self, uid, stype, text, start, stop):
        Position.__init__(self, text, start, stop)

        self.stype = stype.lower()
        assert self.stype in ["chemical", "disease"]

        uids = [Ontology_ID(v) for v in uid.split("|")]
        self.update_uid(frozenset(uids))

    def update_uid(self, new_uid):
        """Update the set of Ontology IDs for this annotation."""
        assert isinstance(new_uid, frozenset), "New UID {0} is not a frozenset!".format(new_uid)
        self.uid = new_uid
        self.has_mesh = self.has_mesh_id()
        self.flat_repr = "|".join([v.flat_repr for v in self.uid])

    def __repr__(self):
        return "<{0}>: '{1}'({2}:{3}) {4}-{5}".format(self.__class__.__name__,
            self.text, self.stype, self.uid, self.start, self.stop)

    def __lt__(self, other):
        """Sort annotations by starting position."""
        if isinstance(other, self.__class__):
            return self.start < other.start

        return NotImplemented

    def has_mesh_id(self):
        """Check if any of the Ontology_IDs are MeSH IDs.

        The gold standard only contains relations between MeSH identified
        concepts.

        E.g., PMID 8638876 (development set)
            The annotation used for a relation is "D002544|-1". No other
            annotations contain D002544 as an identifier.
        """
        for identifier in self.uid:
            if identifier.uid_type == "MESH":
                return True

        return False

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
        self.pmid = int(pmid)
        self.uid = "{0}_{1}".format(pmid, idx)

        # start = position of sentence in the abstract
        Position.__init__(self, text, start, stop)

        # a list of the concept annotations within this sentence
        self.annotations = annotations

        # generate the list of CID and non-CID relations bound to this sentence
        self.poss_relations = self.classify_relations()

    def __repr__(self):
        return "<{0}>: PMID: {1} '{2}'({3}-{4})\nAnnotations: {5}\n".format(
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
            if annot_A.stype == "chemical" and annot_A.has_mesh:
                for annot_B in self.annotations:
                    if annot_B.stype == "disease" and annot_B.has_mesh:
                        # relations between pairs of frozensets
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

class Simple_Rel:
    """A single simple chemical-induced disease relationship.

    Links exactly one chemical Ontology_ID with one disease Ontology_ID. Assumes
    that the chemical and disease are MeSH only.
    """
    def __init__(self, pmid, chemical, disease):
        self.pmid = int(pmid)

        assert isinstance(chemical, Ontology_ID), "Not Ontology_ID"
        assert isinstance(disease, Ontology_ID), "Not Ontology_ID"

        assert chemical.uid_type == "MESH", "Not MeSH ID"
        assert disease.uid_type == "MESH", "Not MeSH ID"

        self.chemical = chemical
        self.disease = disease

    def __repr__(self):
        return "<{0}>: #{1} {2}&{3}".format(self.__class__.__name__,
            self.pmid, self.chemical.flat_repr, self.disease.flat_repr)

class Relation(Base):
    """A single chemical-induced disease relationship.

    Used to compare identifier pairs against the gold standard. Create relations
    from strings with the following format:

    D123456 or MESH:D123456 or MESH:D123456|MESH:D123456

    Concepts linked by a relation are represented as frozensets of Ontology_ID
    objects. Two relations are considered equal if both the chemical and disease
    identifier sets share at least one Ontology_ID in common.
    """
    def __init__(self, pmid, chemical_id, disease_id):
        self.pmid = int(pmid)

        chem_ids = [Ontology_ID(v) for v in chemical_id.split("|")]
        dise_ids = [Ontology_ID(v) for v in disease_id.split("|")]

        chem_ids = [v.flat_repr for v in chem_ids if v.uid_type == "MESH"]
        dise_ids = [v.flat_repr for v in dise_ids if v.uid_type == "MESH"]

        self.uid = "PMID {0}: {1}-{2}".format(pmid, "|".join(chem_ids), "|".join(dise_ids))

        self.chemical_id = frozenset(chem_ids)
        self.disease_id = frozenset(dise_ids)

    def __eq__(self, other):
        """Two relations are considered equal if both the chemical and disease
        identifier sets share at least one Ontology_ID in common.

        This unintuitive definition is used because the gold standard relations
        link two single MeSH identifiers together. However, each identifier may
        only exist as one of many identifiers for a particular concept
        annotation. Therefore it is impossible to determine the origin of a
        relation's MeSH identifier.

        WARNING: This definition of equality breaks transitivity. If there are
        three Relation objects A, B, and C, then if A == B and B == C, then
        A == C is not necessarily true!
        """
        if isinstance(other, self.__class__):
            return (len(self.chemical_id & other.chemical_id) > 0
                and len(self.disease_id & other.disease_id) > 0)

        return NotImplemented

    def __repr__(self):
        return "<{0}>: {1}".format(self.__class__.__name__, self.uid)

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

        self.gold_relations = gold_relations # may be empty when not parsing gold

        self.chemicals, self.diseases = self.get_unique_concepts()

        # split sentences and generate sentence-bound relations
        self.sentences = self.split_sentences()
        self.poss_relations = self.classify_relations()

    def __repr__(self):
        return ("<{0}>: PMID {1}. {2} annotations, {3} gold relations\n"
            "{4} unique chemical ids, {5} unique disease ids\n"
            "{6} sentences".format(self.__class__.__name__,
            self.pmid, len(self.annotations), len(self.gold_relations),
            len(self.chemicals), len(self.diseases), len(self.sentences))
        )

    def has_correct_annotations(self):
        """Check that the annotation indicies produce the correct text snippet.
        """
        full_text = "{0} {1}".format(self.title, self.abstract)
        for annotation in self.annotations:
            fragment = full_text[annotation.start : annotation.stop]
            assert fragment == annotation.text, ("PMID {0} {1} does not match"
                "the text.".format(self.pmid, annotation))

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
        full_text = "{0} {1}".format(self.title, self.abstract)

        # if an abbreviation is included in parentheses, then it should
        # follow the definition annotation immediately
        for i, definition in enumerate(self.annotations[ : -1]):
            if not used[i] and definition.has_mesh:
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
                            and not annot.has_mesh
                            and annot.text == acronym.text):

                            self.annotations[i + 1 + j].update_uid(definition.uid)
                            used[i + 1 + j] = True

    def get_unique_concepts(self):
        """Determine the unique chemical and disease identifiers for this paper.

        Ignores any annotations which are not assigned MeSH identifiers.
        """
        res = defaultdict(set)
        for annotation in self.annotations:
            if annotation.has_mesh:
                res[annotation.stype].add(annotation.uid)

        return (res["chemical"], res["disease"])

    def get_all_possible_relations(self):
        """Return all possible unique chemical-disease identifier pairs."""
        return {(chemical_id, disease_id) for chemical_id in self.chemicals for disease_id in self.diseases}

    def split_sentences(self):
        """Split the abstract into individual sentences, and determine which
        concept annotations reside within each sentence.

        Time complexity: O(N + M) where N is the number of sentences and M is
        the number of annotations.
        """
        all_sentences = [self.title] + split_abstract(self.abstract)

        full_text = "{0} {1}".format(self.title, self.abstract)

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
        sentence_non_cid_rels = set()
        for sentence in self.sentences:
            cid_rels |= sentence.poss_relations[True]
            sentence_non_cid_rels |= sentence.poss_relations[False]

        sentence_non_cid_rels -= cid_rels

        not_sent_bound_rels = all_rels - cid_rels - sentence_non_cid_rels

        assert cid_rels.isdisjoint(sentence_non_cid_rels)
        assert cid_rels.isdisjoint(not_sent_bound_rels)
        assert not_sent_bound_rels.isdisjoint(sentence_non_cid_rels)

        assert (len(self.chemicals) * len(self.diseases)
            == len(cid_rels | sentence_non_cid_rels | not_sent_bound_rels))

        poss_relations = {
            "CID": cid_rels,
            "sentence_non_CID": sentence_non_cid_rels,
            "not_sentence_bound": not_sent_bound_rels
        }
        return poss_relations

    def has_relation(self, potential_relation):
        """Check if the provided relation matches a gold standard relation.

        It is not possible to use a set to check if the relation exists in the
        gold standard because there is no way to make the hashes of two
        Relation objects equal.
        """
        return potential_relation in self.gold_relations

def parse_input(loc, fname, is_gold = True, return_format = "list",
    fix_acronyms = True):
    """Parse a PubTator formatted file and return a set of Paper objects."""
    assert return_format in ["list", "dict"]

    papers = []
    counter = 0
    annotations = []
    relations = []
    for i, line in enumerate(read_file(fname, loc)):
        if len(line) == 0:
            # time to create the paper object
            if is_gold:
                papers.append(Paper(pmid, title, abstract, annotations,
                    relations, fix_acronyms = False))
            else:
                papers.append(Paper(pmid, title, abstract, annotations,
                    [], fix_acronyms = fix_acronyms))

            counter = 0
            annotations = []
            relations = []
        else:
            if 0 <= counter <= 1:
                vals = line.split('|')
                assert len(vals) == 3, "Bad format for line {0}".format(i+1)
            else:
                vals = line.split('\t')

            if counter == 0:
                assert vals[1] == "t", "Bad format for line {0}".format(i+1)
                pmid = int(vals[0])
                title = vals[2]
            elif counter == 1:
                assert vals[1] == "a", "Bad format for line {0}".format(i+1)
                assert int(vals[0]) == pmid
                abstract = vals[2]
            elif is_gold and len(vals) == 4:
                assert int(vals[0]) == pmid and vals[1] == "CID"
                relations.append(Relation(pmid, vals[2], vals[3]))
            else:
                # an annotation
                if len(vals) == 5: # no identifier was assigned
                    vals.append("-1")

                assert 6 <= len(vals) <= 7, "Error on line {0}".format(i+1)
                annotations.append(Annotation(vals[5], vals[4], vals[3], vals[1], vals[2]))

            counter += 1

    if return_format == "dict":
        return {paper.pmid : paper for paper in papers}

    return papers
