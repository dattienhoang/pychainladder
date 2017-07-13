# pychainladder
some basic actuarial methods in python. currently I have some static methods down:
- basic converters between incremental and cumulative triangle types
- Chain Ladder (basic)
- Clark Cape Cod Method (LDF)

in development are some Bayesian stochastic methods from literature:
- Bayesian MCMC for Incurred Loss Triangles (Meyers 2015, CAS Monograph 1)
  - Correlated Chain Ladder
  - Level Chain Ladder
- Bayesian MCMC Models for Incremental Paid Loss Triangles (Meyers 2015, CAS Monograph 1)
  - Correlated Incremental Trend
  - Level Incremental Trend
- Bayesian MCMC Models for Cumulative Paid Loss Triangles (Zhang et al. 2012, J. Royal. Stat. Soc. A)
  - 'non-linear model' (Bayesian MCMC of Clark Cape Cod Model)

References:
[Clark Cape Cod](http://www.casualfellow.com/study-guide.html)
[Bayesian MCMC Reserving Methods](http://www.casact.org/pubs/monographs/index.cfm?fa=meyers-monograph01)
[A Bayesian non-linear model for forecasting insurance loss payments](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-985X.2011.01002.x/full)
