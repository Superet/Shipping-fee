Read me

This folder includes simulation scripts for shipping fee model 

“cond_dens.r”: tests conditional density estimation, compares mixture normal density and nonparametric density. 

“plot_utility.r”: plot utility function for illustration. 

“sim_2constraint.R”: simulate an allocation problem with 2 constraints. 

“sim_prob1.r”: simulate and estimate problem 1 (no shipping fee)

“sim_prob3_0param.r”: problem 3 models total expenditure under shipping fee and search cost. This script explores the roles of parameters in generating to have identification intuition. 

“sim_prob3_1NoCovariates.r”: this script tests the efficient simulation method and tests estimation when restricting the response to shipping fee being constant 1, i.e., with search cost only. 

“sim_prob3_2NoCovariates_2par.r”: tests the identification of response to shipping fee (alpha1) and search cost(alpha2). 

“sim_prob3_3NoCovariates_het.r”: adds heterogeneity to alpha1 and alpha2. 

“sim_prob3_4Covariates_mixture.r”: tests the estimation with covariates, using conditional mixture normal density estimation as the approach in simulated MLE. 

“sim_prob3_4Covariates_np.r”: tests the estimation with covariates, using conditional nonparametric density estimation as the approach in simulated MLE. 

