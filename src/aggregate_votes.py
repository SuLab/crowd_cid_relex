# last updated 2015-06-24 tong shu li

# this is the function for aggregating the votes for crowdflower
# jobs 746297 and 746647
import pandas as pd
from collections import defaultdict

def aggregate_votes(uniq_id, data_frame):
    """
    Given a data frame representing all the unique votes
    for one work unit, aggregates the votes for each of the
    possible choices.

    Returns an unsorted data frame containing the relationships
    with normalized scores.
    """
    rel_id = dict()

    # first map the ids: choice # -> id_pair
    for i in range(5):
        colname = "choice_{0}_ids".format(i)
        assert len(data_frame[colname].unique()) == 1
        rel_id["choice_{0}".format(i)] = data_frame.iloc[0][colname]

    scores = defaultdict(float)
    # increment each relationship pair id by the worker's trust score
    for idx, row in data_frame.iterrows():
        # check that none of the above does not conflict with the other choices
        user_choices = row["chemical_disease_relationships"].split('\n')
        if "none_are_true" in user_choices:
            assert len(user_choices) == 1, idx
            # vote against all other choices
            for i in range(5):
                scores[rel_id["choice_{0}".format(i)]] -= row["_trust"]
        else:
            for choice in user_choices:
                scores[rel_id[choice]] += row["_trust"]

    total_trust = sum(data_frame["_trust"])

    # normalize choices and remove those below zero or empty
    temp = defaultdict(list)
    for id_pair, score in scores.items():
        score /= total_trust
        if score > 0 and id_pair != "empty":
            temp["id_pair"].append(id_pair)
            temp["normalized_score"].append(score)

    df = pd.DataFrame(temp)

    df["uniq_id"] = uniq_id
    assert len(data_frame["_unit_id"].unique()) == 1
    df["unit_id"] = data_frame["_unit_id"].iloc[0]

    return df
