# Tong Shu Li
# Last updated: 2015-11-25

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

    # check to see that class is a binary classifier
    # if there is only one value for the classifier (either 0 or 1),
    # then the ROCR package will throw an error
    # therefore we return either 0 or 1 depending on what the value is
    if (length(unique(class)) == 1)
        return (class[1])

    pred <- prediction(scores, class)

    # draws a ROC curve only if the curve_filename is specified
    if (curve_fname != "")
    {
        perf <- performance(pred, measure = "tpr", x.measure = "fpr")
        png(filename = curve_fname, height = 800, width = 800, bg = "white")

        plot(perf, main = title, type = "o", colorize = TRUE,
            colorkey = TRUE, colorkey.relwidth = 0.5,
            col = rainbow(10))
        abline(a = "0", b = "1", lty = 2)
        dev.off()
    }

    perf <- performance(pred, "auc")
    auc <- slot(perf, "y.values")[[1]]

    return(auc)
}
