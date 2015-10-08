# Tong Shu Li
# Last updated: 2015-10-08

from .F_score import F_score

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
