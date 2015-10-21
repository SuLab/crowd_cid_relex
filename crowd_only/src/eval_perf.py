# Tong Shu Li
# Last updated: 2015-10-21

from collections import defaultdict
from itertools import groupby
import pandas as pd

from .F_score import F_score
from .mesh_filter import filter_relations
from .data_model import Simple_Rel
from .data_model import Ontology_ID

def get_triples(dataframe):
    TRIPLE = ["pmid", "chemical_id", "disease_id"]
    return set(
        dataframe[TRIPLE].apply(
            lambda row: (int(row["pmid"]), row["chemical_id"], row["disease_id"]), axis = 1
        )
    )

def performance(gold, predict, human_readable = False):
    """Calculates precision, recall, and F1 score.

    Given two sets of data objects, calculates as follows:

    A = gold, P = predict

    tp=|A∩P|
    fp=|P-A|
    fn=|A-P|
    The precision (p), recall (r) and f-score (f) are calculated using the typical definitions:

    p=tp/(tp+fp)
    r=tp/(tp+fn)
    f=(2∙p∙r)/(p+r)
    """
    tp = gold & predict
    fp = predict - gold
    fn = gold - predict

    precision = len(tp) / (len(tp) + len(fp))
    recall = len(tp) / (len(tp) + len(fn))
    f1 = F_score(precision, recall)

    if not human_readable:
        return (precision, recall, f1)

    print("# True pos: {0}".format(len(tp)))
    print("# False pos: {0}".format(len(fp)))
    print("# False neg: {0}".format(len(fn)))

    print("Precision: {0}\nRecall: {1}\nF-score: {2}".format(precision, recall, f1))


def official_F_score(score_column, gold_rel_set, dataframe, apply_mesh_filter = False):
    def apply_filter(predict):
        """Apply MeSH Ontology filter."""
        res = set()

        # group by pmid
        temp = sorted(list(predict), key = lambda val: val[0])
        for pmid, group in groupby(temp, lambda val: val[0]):
            rels = [Simple_Rel(pmid, Ontology_ID(chem), Ontology_ID(dise)) for info, chem, dise in group]
            filtered = filter_relations(rels)
            filtered = [(rel.pmid, rel.chemical.flat_repr, rel.disease.flat_repr) for rel in filtered]

            res |= set(filtered)

        return res

    EPSILON = 0.0000001

    res = defaultdict(list)
    for threshold in dataframe[score_column].unique():
        sub = dataframe.query("{0} > {1} or -{2} <= {0} - {1} <= {2}".format(score_column, threshold, EPSILON))

        # grab the relation ids we guessed
        predict = get_triples(sub)

        if apply_mesh_filter:
            predict = apply_filter(predict)

        precision, recall, f1 = performance(gold_rel_set, predict)

        res["recall"].append(recall)
        res["precision"].append(precision)
        res["threshold"].append(threshold)
        res["F_score"].append(f1)

    return pd.DataFrame(res).sort("threshold").reset_index(drop = True)
