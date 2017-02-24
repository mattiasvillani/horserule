# R package for the HorseRule model

The model is described in the paper:

*Tree Ensembles with Rule Structured Horseshoe Regularization*

by Malte Nalenz and <a href="http://mattiasvillani.com">Mattias Villani</a>

*Abstract*: We propose a new Bayesian model for flexible nonlinear regression and classification using tree ensembles. The model is based on the RuleFit approach in Friedman and Popescu (2008) where rules from decision trees and linear terms are used in a L1-regularized regression. We modify RuleFit by replacing the L1-regularization by a horseshoe prior, which is well known to give aggressive shrinkage of noise predictor while leaving the important signal essentially untouched. This is especially important when a large number of rules are used as predictors as many of them only contribute noise. Our horseshoe prior has an additional hierarchical layer that applies more shrinkage a priori to rules with a large number of splits, and to rules that are only satisfied by a few observations. The aggressive noise shrinkage of our prior also makes it possible to complement the rules from boosting in Friedman and Popescu (2008) with an additional set of trees from random forest, which brings a desirable diversity to the ensemble. We sample from the posterior distribution using a very efficient and easily implemented Gibbs sampler. The new model is shown to outperform state-of-the-art methods like RuleFit, BART and random forest on 16 datasets. The model and its interpretation is demonstrated on the well known Boston housing data, and on gene expression data for cancer classification. The posterior sampling, prediction and graphical tools for interpreting the model results are implemented in a publicly available R package.

The paper is available on arXiv: https://arxiv.org/abs/1702.05008

The datasets used in the paper (which are not supplied by R) is included in the folder datasets.

# Example code

Here is an example of how to use the basics of the horserule package to analyze the famous Boston housing data.

    # Loading the library, reading the data
    library(horserule)
    data(Boston)

    # Training/Testing split of the data
    N = nrow(Boston)
    train = sample(1:N, 500)
    Xtrain = Boston[train,-14]
    ytrain = Boston[train, 14]
    Xtest = Boston[-train, -14]
    ytest = Boston[-train, 14]

    # Main function call (variable scaling performed internally)
    hrres = HorseRuleFit(X=Xtrain, y=ytrain,
		# MCMC settings
		thin=1, niter=1000, burnin=100,
		# Parameters for the rule generation process
		L=5, S=6, ensemble = "both", mix=0.3, ntree=1000,
		# Model parameters. Data is scaled so no intercept needed.
		intercept=F, linterms=lin, ytransform = "log",
		# Hyperparameters for the rule structured prior
		alpha=1, beta=2, linp = 1, restricted = 0)

    # Check model performance by predicting holdout cases
    pred = predict(hrres, Xtest)
    sqrt(mean((pred-ytest)^2)))

    # Find most important rules
    importance_hs(hrres)

	# Compute variable importance
	variable_importance(hrres)

	# Monitor the complexity of the rules
	complexity_plot(hrres)
