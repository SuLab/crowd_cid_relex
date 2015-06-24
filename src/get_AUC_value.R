# last updated 2015-06-24 tong shu li

# calculates the AUC value for a ROC curve
# without actually generating a plot of the curve
# since rpy2 cannot convert pandas dataframes
# to R dataframes, the data frame has to be read
# from a text file. The score_col and class_col
# variables determine the column names to be used
# for generating the ROC curve

suppressMessages(library(ROCR))

AUC_value <- function(fname, score_col, class_col)
{
    df <- read.table(fname, sep = '|', header = TRUE)

    scores <- df[ , score_col]
    class <- df[ , class_col]

    pred <- prediction(scores, class)

    temp <- performance(pred, "auc")
    auc <- slot(temp, "y.values")[[1]]

    return(auc)
}
