"""
Tong Shu Li
Created on Tuesday, 2015-07-28
Last updated 2015-07-28

This file contains the functions for determining the F-score,
precision, and recall for a general dataframe.
"""
from collections import defaultdict
import pandas as pd

def F_score(precision, recall):
    return 2 * (precision * recall) / (precision + recall)

def all_F_scores(score_column, class_column, data_frame):
    """
    Given a dataframe of results, calculates the F-score,
    precision, and recall for all possible unique
    threshold levels of the score column.

    score_column = name of column containing the threshold rating score
    class_column = name of column containing 0 if not true according to
        the gold standard, and 1 if true
    """
    EPSILON = 0.0000001

    res = defaultdict(list)

    for threshold in data_frame[score_column].unique():
        sub = data_frame.query("{0} > {1} or -{2} <= {0} - {1} <= {2}".format(
            score_column, threshold, EPSILON)
        )

        recall = sum(sub[class_column]) / sum(data_frame[class_column])
        precision = sum(sub[class_column]) / len(sub)

        f_score = F_score(precision, recall)

        res["recall"].append(recall)
        res["precision"].append(precision)
        res["F_score"].append(f_score)
        res["threshold"].append(threshold)

    return pd.DataFrame(res)

def max_F_score(score_column, class_column, data_frame):
    """
    Finds the maximum F-score for a given dataframe with
    results. Also returns the precision, recall, and
    threshold at which the maximum F-score is achieved.
    """
    res = all_F_scores(score_column, class_column, data_frame)
    idx = res["F_score"].idxmax()
    return res.loc[idx]
