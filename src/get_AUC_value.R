# last updated 2015-06-25 tong shu li

# calculates the AUC value for a ROC curve
# if a filename and title are provided,
# then the function also actually plots the curve

# input is taken from a file

# score_col and class_col determine the
# column names used for generating the curve

suppressMessages(library(ROCR))

AUC_value <- function(fname, score_col, class_col,
    curve_fname = "", title = "")
{
    df <- read.table(fname, sep = '|', header = TRUE)

    scores <- df[ , score_col]
    class <- df[ , class_col]

    pred <- prediction(scores, class)

    # draws a ROC curve only if the curve_filename is specified
    if (curve_fname != "")
    {
        perf <- performance(pred, measure = "tpr", x.measure = "fpr")
        png(filename = curve_fname, height = 800, width = 800, bg = "white")

        plot(perf, main = title, type = "o", col = rainbow(10))
        abline(a = "0", b = "1", lty = 2)
        dev.off()
    }

    perf <- performance(pred, "auc")
    auc <- slot(perf, "y.values")[[1]]

    return(auc)
}
