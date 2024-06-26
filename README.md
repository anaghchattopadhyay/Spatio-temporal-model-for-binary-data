# Spatio-temporal-model-for-binary-data

## Code for paper: "A spatio-temporal model for binary data and its application in analyzing the direction of COVID-19 spread"

### About the paper

It is often of primary interest to analyze and forecast the levels of a continuous phenomenon as a categorical variable. In this paper, we propose a new spatio-temporal model to deal with this problem in a binary setting, with an interesting application related to the COVID-19 pandemic, a phenomena that depends on both spatial proximity and temporal auto-correlation. 

Our model is defined through a hierarchical structure for the latent variable, which corresponds to the probit-link function. The mean of the latent variable in the proposed model is designed to capture the trend, the seasonal pattern as well as the lagged effects of relevant regressors. The covariance structure of the model is defined as an additive combination of a zero-mean spatio-temporally correlated process and a white noise process. The parameters associated with the space-time process enable us to analyze the effect of proximity of two points with respect to space or time and its influence on the overall process. 

For estimation and prediction, we adopt a complete Bayesian framework, along with suitable prior specifications and utilize the concepts of Gibbs sampling. Using the county-level data from the state of New York, we show that the proposed methodology provides superior performance than benchmark techniques. We also use our model to devise a novel mechanism for predictive clustering which can be leveraged to develop localized policies.
