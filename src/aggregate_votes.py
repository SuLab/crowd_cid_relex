# last updated 2015-06-27 Tong Shu Li
"""
Aggregates votes for individual work units.
"""
import pandas as pd
from collections import defaultdict

def aggregate_votes(uniq_id, data_frame, penalize_unchosen = False)
    """
    Given a data frame representing all the unique votes
    for one work unit, aggregates the votes for each of the
    possible choices.

    Returns an unsorted data frame containing the relationships
    with normalized scores.

    Two aggregation schemes can be used:
      1. Choices which are unchosen get no change to their score.
      2. Choices which are unchosen get the trust score subtracted.
    """
    scores = defaultdict(float)
    for idx, row in data_frame.iterrows():
        user_choices = row["chemical_disease_relationships"].split('\n')
        if "none_are_true" in user_choices:
            assert len(user_choices) == 1, idx
            for i in range(5):
                scores["choice_{0}".format(i)] -= row["_trust"]
        else:
            for i in range(5):
                choice = "choice_{0}".format(i)
                if choice in user_choices:
                    score[choice] += row["_trust"]
                else:
                    score[choice] += (-1 if penalize_chosen else 0) * row["_trust"]

    rel_id = dict()
    # map the ids: choice # -> id_pair
    for i in range(5):
        colname = "choice_{0}_ids".format(i)
        assert len(data_frame[colname].unique()) == 1
        rel_id["choice_{0}".format(i)] = data_frame.iloc[0][colname]

    total_trust = sum(data_frame["_trust"])

    # normalize scores
    temp = defaultdict(list)
    for choice, score in scores.items():
        score /= total_trust
        if rel_id[choice] != "empty":
            temp["id_pair"].append(rel_id[choice])
            temp["normalized_score"].append(score)

    df = pd.DataFrame(temp)

    df["uniq_id"] = uniq_id
    assert len(data_frame["_unit_id"].unique()) == 1
    df["unit_id"] = data_frame["_unit_id"].iloc[0]

    return df
