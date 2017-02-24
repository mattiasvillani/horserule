# horserule
R package for the HorseRule model.

library(horserule)
data(Boston)

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
