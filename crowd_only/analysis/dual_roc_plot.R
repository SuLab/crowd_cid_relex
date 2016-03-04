# Tong Shu Li
# Last updated: 2015-12-07
# Plot the dual ROC plot used for the Database BioCreative paper

library(ROCR)

# no NER error filtering
df1 <- read.table("true_roc.tsv", sep = '\t', header = T)
# with the NER error filter applied
df2 <- read.table("no_ner_roc.tsv", sep = '\t', header = T)

scores <- df1[ , "num_votes"]
class <- df1[ , "in_gold"]
pred1 <- prediction(scores, class)

scores <- df2[ , "num_votes"]
class <- df2[ , "in_gold"]
pred2 <- prediction(scores, class)

perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
perf2 <- performance(pred2, measure = "tpr", x.measure = "fpr")

png(filename = "crowd_testset_perf_roc.png",
    units = "in", res = 600, bg = "white",
    height = 10, width = 10)

par(cex.axis = 1.5, cex.lab = 2,
    mar = c(5, 5, 1, 3),
    pty = "s")

plot(perf1, type = "o", lwd = 5, pch = 4)

plot(perf2, add = TRUE, type = "o", lwd = 5, lty = 3)

abline(a = "0", b = "1", lty = 2)

legend("bottomright", c("Without NER filter", "With NER filter"),
       lty = c(1, 3), lwd = 5, cex = 2)

dev.off()
