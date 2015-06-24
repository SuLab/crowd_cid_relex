# last updated 2015-06-24 tong shu li

# This is the Python wrapper for the R get_AUC_value.R
# function. Since pandas dataframes can't be converted
# to R dataframes with rpy2, we have to unfortunately
# write the dataframe to a text file first

import pandas as pd
import rpy2.robjects as robjects
import tempfile

def get_AUC_value(data_frame, score_col, class_col):
    """
    Given a Pandas dataframe with the names of the
    score and class columns, calls the R function for
    determining the AUC value of the ROC curve.

    Methods:
      1. Makes a temporary file (default is in /tmp)
      2. Writes the pandas dataframe to the temporary file
      3. Calls R function using rpy2
      4. R function reads from temporary file and calculates AUC
      5. Temporary file is deleted
      6. AUC value is returned to the calling Python function
    """
    with tempfile.NamedTemporaryFile(mode = "w", delete = True) as fout:
        fname = fout.name

        data_frame.to_csv(fname, sep = "|", index = False)

        R = robjects.r
        R.source("get_AUC_value.R")

        auc = R.AUC_value(fname, score_col, class_col)

    return auc[0]
