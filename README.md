## SG GxE Prediction

Pipeline and helper functions used by the LSU Small Grains breeding program to analyze historic SunGrains data sets using environmental data for predicting merit in groups of counties/parishes. 

The associated files for pulling in data, running model comparisons/cross validation, and running predictions have been iterated on but are a little disorganized and messy. Here I'll attempt to clean them up a bit and make them more repeatable! (With an eye to making the process of changing the inputs and genotypes to predict easier from year-to-year in the future).

Original modeling framework for the BGLR models incorporating interaction matrices between genotype PCs/data and environmental PCs or data from Rogers *et al.* 2022 (https://academic.oup.com/g3journal/article/12/2/jkab440/6486423).

