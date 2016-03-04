"""
Tong Shu Li
Last updated 2015-08-11
"""
import pandas as pd
from collections import defaultdict
import pycountry

def get_country_name(code):
    return pycountry.countries.get(alpha3 = code).name.split(",")[0]

def build_worker_profiles(raw_data):
    """
    Build a profile for every worker that tried the task.
    """
    res = defaultdict(list)

    metadata = ["_trust", "_ip", "_channel"]
    ans_choices = ["yes_direct", "yes_indirect", "no_relation", "ner_mistake"]

    for worker_id, group in raw_data.groupby("_worker_id"):
        test_resp = group.query("_golden")
        work_resp = group.query("~_golden")

        res["worker_id"].append(worker_id)
        res["test_ques_seen"].append(len(test_resp["uniq_id"].unique()))
        res["work_units_seen"].append(len(work_resp["uniq_id"].unique()))

        res["country"].append(get_country_name(test_resp["_country"].iloc[0]))

        for metadata_col in metadata:
            res[metadata_col.lstrip("_")].append(test_resp[metadata_col].iloc[0])

        for work_type, resp_data in zip(["test", "work"], [test_resp, work_resp]):
            time_series = determine_time_taken(resp_data) # time per page
            time_series /= WORK_UNITS_PER_PAGE # time per work unit

            stats = time_stats(time_series)
            for i, name in enumerate(["min", "median", "max"]):
                res["{0}_{1}_time_per_unit".format(work_type, name)].append(stats[i])

            # look at the response distributions
            for ans_choice in ans_choices:
                temp = resp_data.query("verify_relationship == '{0}'".format(ans_choice))
                res["{0}_{1}".format(work_type, ans_choice)].append(len(temp["uniq_id"].unique()))

    return pd.DataFrame(res)
