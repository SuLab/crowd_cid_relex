# last updated 2015-06-04 toby
import pandas as pd
import os

def filter_data(settings):
    """
    Filters raw CrowdFlower output down to data we care about.
    Formatting the data is left to another program.
    """
    data = pd.read_csv(os.path.join(settings["loc"], settings["fname"]), sep = ",")

    if settings["data_subset"] == "gold":
        data = data.query("_golden")
    elif settings["data_subset"] == "normal":
        data = data.query("~_golden")

    data = (data.query("{0} <= _trust <= {1}".
        format(settings["min_accuracy"], settings["max_accuracy"])))

    return data
