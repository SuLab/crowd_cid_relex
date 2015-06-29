# last updated 2015-06-29 tong shu li
"""
This is the Python wrapper for get_AUC_value.R
The R function always returns the AUC value, and
only draws the curve if the filename and title
are specified.
"""
import os
import pandas as pd
import rpy2.robjects as robjects
import tempfile

def get_AUC_value(data_frame, score_col, class_col,
    curve_fname = "", plot_title = ""):
    """
    Given a Pandas dataframe with the names of the
    score and class columns, calls the R function for
    determining the AUC value of the ROC curve.

    Methods:
      1. Makes a temporary file (default is in /tmp)
      2. Writes the pandas dataframe to the temporary file
      3. Calls R function using rpy2
      4. R function reads from temporary file and calculates AUC,
        and optionally draws the ROC curve.
      5. Temporary file is deleted
      6. AUC value is returned to the calling Python function
    """
    with tempfile.NamedTemporaryFile(mode = "w", delete = True) as fout:
        fname = fout.name

        data_frame.to_csv(fname, sep = "|", index = False)

        R = robjects.r
        R.source(os.path.realpath(__file__).replace(".py", ".R"))
        auc = R.AUC_value(fname, score_col, class_col, curve_fname, plot_title)

    return auc[0]
