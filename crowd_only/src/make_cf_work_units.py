# Tong Shu Li
# Created on: 2015-12-21
# Last updated: 2015-12-21
"""
Create CrowdFlower work units.
"""
from collections import defaultdict
import pandas as pd

from .make_sections import create_sections

def add_simple_tag(tag_name, tag_class, text):
    """Put a HTML tag around a snippet of text."""
    return "<{0} class=\"{1}\">{2}</{0}>".format(tag_name, tag_class, text)

def highlight_concepts(text, breaks):
    """
    Inserts HTML tags around the pieces of text
    which need to be highlighted in a string.
    """
    breaks = sorted(breaks, key = lambda x: x[0])

    final = []
    for i in range(len(breaks) - 1):
        s = text[breaks[i][0] : breaks[i+1][0]]
        if breaks[i][1] != "n":
            s = add_simple_tag("span", breaks[i][1], s)

        final.append(s)

    return "".join(final)

def highlight_text(text, offset, uniq_spans):
    """
    Given a string and the annotations which fall
    within this string, highlights the concepts.
    """
    # index of break, type of break (n = nothing)
    breaks = [(0, "n"), (len(text), "n")]

    for span in uniq_spans:
        breaks.append((span.start - offset, span.stype))
        breaks.append((span.stop - offset, "n"))

    return highlight_concepts(text, breaks)

def grab_names(annotations):
    """Determines the unique names of the annotations."""
    used_names = defaultdict(set) # lower case set of used names (to avoid repeats)
    real_name = defaultdict(set) # set of unique names verbatim (to preseve capitalization)
    for annotation in annotations:
        if annotation.text.lower() not in used_names[annotation.stype]:
            used_names[annotation.stype].add(annotation.text.lower())
            real_name[annotation.stype].add(annotation.text)

    return real_name

def process_sentence(chem, dise, sentence):
    spans = [annot for annot in sentence.annotations if not {chem, dise}.isdisjoint(annot.uid.uid)]

    form_text = highlight_text(sentence.text, sentence.start, spans)

    names = grab_names(spans)

    if len(set(v.stype for v in spans)) == 2:
        form_text = add_simple_tag("span", "sentence", form_text)

    return (names, form_text)

def create_work_units(dataset):

    def get_rels(paper):
        for origin, rels in paper.poss_relations.items():
            if origin != "CID":
                for rel in rels:
                    yield (rel, origin)

    cid = dict()

    res = defaultdict(list)
    for pmid, paper in dataset.items():
        cid[pmid] = paper.poss_relations["CID"]

        for (chem, dise), origin in get_rels(paper):

            # individually highlight each sentence and grab names
            sentences = []
            names = defaultdict(set)
            for s in paper.sentences:
                snippets, form_sent = process_sentence(chem, dise, s)
                sentences.append(form_sent)
                for stype in ["chemical", "disease"]:
                    names[stype] |= snippets[stype]

            # join together

            form_title = sentences[0]
            form_body = " ".join(sentences[1:])

            form_body = create_sections(form_body)

            res["pmid"].append(pmid)
            res["chemical_id"].append(chem.flat_repr)
            res["disease_id"].append(dise.flat_repr)
            res["rel_origin"].append(origin)

            res["chemical_name"].append(add_simple_tag("span", "chemical",
                "/".join(names["chemical"])))

            res["disease_name"].append(add_simple_tag("span", "disease",
                "/".join(names["disease"])))

            res["form_title"].append(form_title)
            res["form_body"].append(form_body)

    res = pd.DataFrame(res)
    res["uniq_id"] = pd.Series(["refine_try_1_{}".format(i) for i in res.index])

    return (cid, res)
