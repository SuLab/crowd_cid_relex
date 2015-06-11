Readme. 

We already have the results of your proposal to measure the impact of training set size on performance. We have blocked out the 500 development articles into groups of 100, 200, 300, 400 and 500, to train those 5 models and test each on the 500 abstracts in the training set, and vice versa. I have attached an image ("BC5_eval.png") representing the different F-scores obtained in each evaluation.  We can see from the graphs that the performance improves when the number of examples for training is higher.

In addition, a set of different evaluations were carried out with the EUADR-crowd and BC5 TraningSet. The  "evaluations.jpeg" file shows the results of the evaluations:

    First row: the 10-fold Cross-Validation (10F-CV) of the EUADR-crowd corpus.
    Second row: the results obtained from the model trained with the EUADR-crowd corpus against the BC5 TraningSet.
    Third row: similar to the first row, but taking as true association only the sentences annotated as "may_cause".
    Fourth row: similar to the second row, but taking as true association only the sentences annotated as "may_cause".
    Fifth row: the 10F-CV of the BC5 TraningSet.

We think that the number of examples from the EUADR-crowd corpus is not high to measure the impact of training set size on performance.