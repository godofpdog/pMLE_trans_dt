# The pseudo maximum likelihood estimators of semiparametric transformation models with bouby truncated models

## douby truncated data
* data with both left-truncation and right-truncation.

### truncation
* censored : partial information about the value of a variable (only know it is beyond some boundary).
* truncation : data set does not include observations in the analysis that are beyond a boundary value.
* We can treat the truncation as bias sampling. 
* If we do not consider the truncation, it will be overfitting for the population we really interest in.
### left truncation
If the recuitment of patients continues after the onset time of a stuty, patients who have aready experienced the event are often excluded from the study, which results the left trunction of the event time we are insterest in.
### right trunction
* retrospective sampling criterion
* long latent period

## semiparametric transformation models
* Cox proportion hazard model is the most popular semiparametric model in survival analysis, but the proportion hazard assumption can be violated in some situations.
* We can treat the semiparametric transformation model as a generalized Cox ph model.
* In application, we just find a good transformation model to fit the data by some metrics just like tuning hyper-parameters in machine learning area.

## the API
* simulation 
* application

## simulation : 
Simulation is the way to find out the performance of entire algorithm on different situations. Generating data by a known parameters of model, we can get some metrics like 'rmse' of the estimation of weights in the model                         (see '/simulation/output/summary.txt'). That is, we controll the situations by setting parameters then evaluate the algorithm's performane and decide whether to use this method.

### *RUN.R* : run the algorithm by setup some params.
* n : sample size of the simulation
* REP : replication of whole simulation.
* d0 : parameter to control the truncation rate.
* sbeta : initial values of beta related to the event time distribution.
* rho : the parameter of transformation model. (log trans, BoxCox trans)
* R.in : inverse function of baseline hazerd.
* md : model of the gap time.
* Gnt_est : method to estimating the truncated distibution.
### *continue.R* : continue the programs if it shutdown.
### *setup.R* help us to find the params to generate the data of wanted situation.
*  see 'simulation/programs/data_gen.R' for details of data generating.
