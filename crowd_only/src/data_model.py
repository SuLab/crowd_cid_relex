# Tong Shu Li
# First written: 2015-07-02
# Last updated: 2015-08-17
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
from itertools import islice

from lingpipe.file_util import read_file
from lingpipe.split_sentences import split_abstract

def is_MeSH_id(uid):
    return len(uid) == 7 and uid[0] in ["C", "D"]

class Ontology_ID:
    """
    A single identifier from an existing ontology.
    """
    def __init__(self, text):
        """
        Given the raw text representation, splits and
        formats accordingly.
        """
        if ":" in text:
            vals = text.split(':')
            assert len(vals) == 2, "ID {0} is misformatted!".format(text)
            self.uid_type = vals[0]
            self.uid = vals[1]

            if self.uid_type == "MESH":
                assert is_MeSH_id(self.uid), "MeSH ID mismatch for {0}".format(text)
        else:
            self.uid = text
            if is_MeSH_id(self.uid):
                self.uid_type = "MESH"
            else:
                self.uid_type = "unknown"

    def flat_repr(self):
        """
        A flat string representation of this ontology ID.
        """
        return "{0}:{1}".format(self.uid_type, self.uid)

    def __repr__(self):
        return "<{0}>: {1}:{2}".format(self.__class__.__name__,
            self.uid_type, self.uid)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.uid_type == other.uid_type and self.uid == other.uid

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash("{0}:{1}".format(self.uid_type, self.uid))

class Annotation:
    """
    A single mention of a concept in a piece of text.
    Annotation positions are indexed to the abstract.
    """
    def __init__(self, uid, stype, text, start, stop):
        uids = map(lambda v: Ontology_ID(v), uid.split("|"))

        self.uid = frozenset(uids)
        self.has_mesh = self.has_mesh_id()

        self.stype = stype.lower()
        assert self.stype in ["chemical", "disease"]

        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop, "Annotation {0} indicies reversed!".format(self.uid)
        assert len(text) == self.stop - self.start, "Annotation {0} length mismatch!".format(self.uid)

    def flat_repr(self):
        """
        A flat representation for printing.
        """
        return "|".join(map(lambda v: v.flat_repr(), self.uid))

    def __repr__(self):
        return "<{0}>: '{1}'({2}:{3}) {4}-{5}".format(self.__class__.__name__,
            self.text, self.stype, self.uid, self.start, self.stop)

    def __cmp__(self, other):
        """
        When sorting Annotation objects, make sure they
        are in ascending order of starting position.
        """
        if hasattr(other, "start"):
            return self.start.__cmp__(other.start)

    def has_mesh_id(self):
        """
        Checks if any of the identifiers used for this annotation are
        MeSH. Only identifiers which have at least one MeSH component
        are included in the gold relations.

        E.g., PMID 8638876 (development)
            The annotation id is "D002544|-1", and two of the gold
            relations have "D002544" as the disease. No other annotation
            contains D002544.
        """
        for identifier in self.uid:
            if identifier.uid_type == "MESH":
                return True

        return False

class Sentence:
    """
    A single sentence from a paper.

    Sentence-bound relationships are generated
    by this object.

    The set of chemical-disease relations that need to
    be made into a work unit is the non-CID relations
    minus the CID relations true at the abstract level.
    """
    def __init__(self, pmid, idx, text, start, stop, annotations):
        self.pmid = int(pmid)
        self.uid = "{0}_{1}".format(pmid, idx)
        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop, "Sentence {0} indicies reversed!".format(self.uid)
        assert len(text) == self.stop - self.start, "Sentence {0} length mismatch!".format(self.uid)

        # a list of the concept annotations within this sentence
        self.annotations = annotations

        # generate the list of CID and non-CID relations bound to this sentence
        self.poss_relations = self.classify_relations()

    def __repr__(self):
        return "<{0}>: PMID: {1} '{2}'({3}-{4})\nAnnotations: {5}\n".format(
            self.__class__.__name__, self.pmid, self.text,
            self.start, self.stop, self.annotations)

    def is_CID_relation(self, chemical, disease):
        """
        Given two Annotations within this sentence,
        determine if the pair of annotations follows the
        CID structure.
        """
        return ((chemical.stop < disease.start)
            and (disease.start - chemical.stop <= 15)
            and ("induce" in self.text[chemical.stop - self.start :
                disease.start - self.start].lower())
        )

    def classify_relations(self):
        """
        This function generates all the unique chemical-disease
        identifier pairs of annotations contained within this
        sentence, and classifies them into two groups:
            1. Those which follow a '[chemical]-induced [disease]' (CID)
                structure.
            2. Those which do not follow the CID structure.

        The two sets of identifier pairs are mutually exclusive.

        Only MeSH identifiers will be generated for CID relation
        verification, since the gold standard will never use any
        other ontology identifier. This has been confirmed with
        Zhiyong. Therefore "-1" identifiers can be skipped, as can
        CHEBI and other ontology identifier types.
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

class Relation:
    """
    A single chemical-induced disease relationship.
    Used to compare identifier pairs against the gold.

    Relation objects should be created from strings with the following
    format:

    D123456 or MESH:D123456 or MESH:D123456|MESH:D123456

    Relation identifiers are represented as frozensets of
    Ontology_ID objects. Two relations are considered equal
    if the chemical and disease identifier set intersection
    is at least one.

    When outputting for the gold, relations should only use one MESH id
    for both chemical and disease, instead of multiple as for the
    annotations.

    All single ontology ids in the frozenset should be MeSH.
    """
    def __init__(self, pmid, chemical_id, disease_id):
        self.pmid = int(pmid)

        chem_ids = map(lambda v: Ontology_ID(v), chemical_id.split("|"))
        dise_ids = map(lambda v: Ontology_ID(v), disease_id.split("|"))

        chem_ids = filter(lambda v: v.uid_type == "MESH", chem_ids)
        dise_ids = filter(lambda v: v.uid_type == "MESH", dise_ids)

        chem_uids = map(lambda v: v.uid, chem_ids)
        dise_uids = map(lambda v: v.uid, dise_ids)

        self.uid = "PMID {0}:{1}-{2}".format(pmid, "|".join(chem_uids), "|".join(dise_uids))

        self.chemical_id = frozenset(chem_ids)
        self.disease_id = frozenset(dise_ids)

    def __eq__(self, other):
        """
        Equal if the intersection of both the chemical
        and disease identifier frozensets is at least
        one element large.

        This is because the gold relations only use
        a pair of single MeSH ids, despite the fact
        that the annotations use complexed MeSH ids.

        WARNING:
            Defining the equals function in this manner
            breaks transitivity. That is, if we have
            three objects A, B, and C, then if
                A == B and B == C, then
                A == C IS NOT TRUE!!
        """
        if isinstance(other, self.__class__):
            return (len(self.chemical_id & other.chemical_id) > 0
                and len(self.disease_id & other.disease_id) > 0)

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "<{0}>: {1}".format(
            self.__class__.__name__, self.uid)

class Paper:
    """
    A single academic publication.
    Contains:
        1. The PubMed identifier.
        2. The title as a string.
        3. The abstract as a string.
        4. A list of all chemical and disease annotations
            in the title and abstract sorted in increasing
            order of starting index.
        5. A potentially empty list of gold standard CID
            relations.
        6. A set of unique chemical identifiers.
        7. A set of unique disease identifiers.
        8. A list of Sentences containing both the title
            and body of the abstract. The first sentence
            is the title. Each sentence contains the
            annotations and relations constrained to that
            particular sentence.
        9. A set of all the potential chemical-disease
            relations grouped into three mutually exclusive
            categories:
            - CID relations
            - Non-CID sentence-bound relations
            - Non-sentence bound relations
            The sum of relations in all three groups should
            equal the number of unique chemical IDs times
            the number of unique disease IDs.
        10. A function for resolving acronyms.
    """
    def __init__(self, pmid, title, abstract, annotations, gold_relations = []):
        self.pmid = int(pmid)
        self.title = title
        self.abstract = abstract

        self.annotations = sorted(annotations)
        assert self.has_correct_annotations()
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
        """
        Checks that the paper's annotations match the
        stated positions in the text.
        """
        text = "{0} {1}".format(self.title, self.abstract)
        for annotation in self.annotations:
            assert text[annotation.start : annotation.stop] == annotation.text, (
                "Annotation {0} in PMID {1} does not match the text.".format(annotation, self.pmid))

        return True

    def resolve_acronyms(self):
        """
        This function tries to resolve acronyms.
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

                            self.annotations[i + 1 + j].uid = definition.uid
                            used[i + 1 + j] = True

    def get_unique_concepts(self):
        """
        Determines the unique identifiers of chemicals and diseases
        belonging to this paper.

        Ignores any annotations which are not MeSH identified.
        """
        res = defaultdict(set)
        for annotation in self.annotations:
            if annotation.has_mesh:
                res[annotation.stype].add(annotation.uid)

        return (res["chemical"], res["disease"])

    def get_all_possible_relations(self):
        """
        Returns all possible unique drug-disease combinations
        as a set for this paper.
        """
        return {(chemical_id, disease_id) for chemical_id in self.chemicals for disease_id in self.diseases}

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
        M = len(self.annotations)
        for i, sentence in enumerate(all_sentences):
            # The sentence splitter isn't perfect. It recognizes "i.v." as a
            # sentence. Since there can be multiple instances of "sentences"
            # like "i.v." (e.g., PMID 10840460), we need to make sure that
            # we are checking for the first instance starting at the current
            # position (since find always finds the first instance otherwise).
            assert full_text.find(sentence, sent_idx) == sent_idx, (
                "PMID {0} sentence '{1}' does not match text!".format(self.pmid, sentence))

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
        """
        Takes all the possible relations for this abstract and
        splits them into three mutually exclusive groups:
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
        return potential_relation in self.gold_relations

def parse_input(loc, fname, is_gold = True, return_format = "list"):
    """
    Reads a given file and returns a list of Paper
    objects.
    """
    assert return_format in ["list", "dict"]
    if return_format == "list":
        papers = []
    else:
        papers = dict()

    counter = 0
    annotations = []
    relations = []
    for i, line in enumerate(read_file(fname, loc)):
        if len(line) == 0:
            # time to create the paper object
            if return_format == "list":
                papers.append(Paper(pmid, title, abstract, annotations, relations))
            else:
                papers[pmid] = Paper(pmid, title, abstract, annotations, relations)

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
                pmid = int(vals[0])
                title = vals[2]
            elif counter == 1:
                assert vals[1] == "a"
                assert int(vals[0]) == pmid
                abstract = vals[2]
            elif is_gold and len(vals) == 4:
                assert vals[1] == "CID"
                relations.append(Relation(pmid, vals[2], vals[3]))
            else:
                assert 5 <= len(vals) <= 7, "Error on line {0}".format(i)

                if len(vals) == 5: # output of tmChem, should be "-1" identifier, but is blank
                    vals.append("-1")

                assert 6 <= len(vals) <= 7, "Error on line {0}".format(i)
                annotations.append(Annotation(vals[5], vals[4], vals[3], vals[1], vals[2]))

            counter += 1

    return papers
