# Tong Shu Li
# Last updated 2015-07-30
"""
Takes raw CrowdFlower data and groups data together by
individual work units. Aggregates the votes based on
a specific column and creates a summary data frame.
"""
import pandas as pd

def aggregate_votes(column_name, data_frame, mapping = None):
    """
    Given all of the human responses for one work unit,
    aggregates the results based on the column you give it.

    For each possible choice, it calculates:
        1. Total number of votes for that choice.
        2. Sum of trust scores of people who chose that choice.

    Returns an unsorted data frame of the results.

    Confidence score is the sum of the trust score of
    individual workers.

    If given a mapping, will convert answers.
    """
    # number of votes and sum trust score for each unique choice
    res = dict()
    for value, group in data_frame.groupby(column_name):
        conf_score = sum(group["_trust"])
        num_votes = len(group) # assuming unique workers for each work unit

        if mapping is not None and value in mapping: # only convert those which need to be changed
            value = mapping[value]

        if value in res:
            res[value] = map(add, res[value], [conf_score, num_votes])
        else:
            res[value] = [conf_score, num_votes]

    temp = [[key] + value for key, value in res.items()] # make rows of data frame
    return pd.DataFrame(temp, columns = [column_name, "conf_score", "num_votes"])

def aggregate_results(id_column, agg_column, data, reporting_method, metadata_columns,
                      reporting_value = None, mapping = None):
    """
    Collect votes together for each individual work unit.
    Choose one of two result reporting methods:
        1. Majority vote
        2. Positive signal only

    reporting_method: {"majority_vote", "positive_signal_only"}
    if reporting_method == "positive_signal_only":
        reporting_value is the label of the positive value in the aggregation column
    """
    assert reporting_method in ["majority_vote", "positive_signal_only"], "Wrong result reporting method!"

    result = []
    for identifier, group in data.groupby(id_column):
        # group the votes together and convert if necessary
        ans = aggregate_votes(agg_column, group, mapping)

        if reporting_method == "positive_signal_only":
            positive_result = ans.query("{0} == '{1}'".format(agg_column, reporting_value))
            if positive_result.empty:
                # no worker chose the positive choice, so create an empty entry
                temp = pd.DataFrame([[reporting_value, 0, 0]], columns = [agg_column, "conf_score", "num_votes"])
                ans = ans.append(temp)

        # figure out what percentage of the workers chose each choice
        total_vote_score = sum(ans["conf_score"])
        ans.loc[:, "percent_agree"] = ans.loc[:, "conf_score"] / total_vote_score

        ans.insert(0, id_column, identifier)

        # copy the metadata columns over (e.g., pmid, unit id)
        for column_name in metadata_columns:
            assert len(group[column_name].unique()) == 1
            value = group[column_name].iloc[0]
            ans.loc[:, column_name.lstrip('_')] = value

        # sort results
        # decreasing confidence score if doing majority vote
        # put positive result at the top if doing positive result method

        if reporting_method == "majority_vote":
            ans = ans.sort(["conf_score"], ascending = False)
            result.append(ans)
        else:
            # put positive result at the top
            positive_result = ans.query("{0} == '{1}'".format(agg_column, reporting_value))
            temp = ans.query("{0} != '{1}'".format(agg_column, reporting_value))
            temp = temp.sort(["conf_score"], ascending = False)

            result.append(pd.concat([positive_result, temp]))

    return pd.concat(result)
