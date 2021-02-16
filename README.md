# Bayesian-Inference-of-Network-Structure-with-Noisy-Data

BAYESIAN STATISTIC PROJECT A.Y. 2020-2021

Our project is about Bayesian inference for network structure with noisy data.  
We make inference about data Xij that represent interaction between students. From this data we want to infer a network, with 0,1 value, that indicates friendship links between the students.
We used four models:

* __Poisson Model - Erdős-Rényi Network Model__
* __Negative Binomial Model - Erdős-Rényi Network Model__
* __Negative Binomial Model - Stochastic Block Model__
* __Poisson Model - Multi-edge Model__

# Poisson Model - Erdős-Rényi Network Model
![Poisson Model - Erdős-Rényi Network Model](https://github.com/NoemiRossi/Bayesian-Inference-of-Network-Structure-with-Noisy-Data/blob/main/Images/p_er_m.png)
# Negative Binomial Model - Erdős-Rényi Network Model
![Negative Binomial Model - Erdős-Rényi Network Model](https://github.com/NoemiRossi/Bayesian-Inference-of-Network-Structure-with-Noisy-Data/blob/main/Images/nb_er_m.png)
# Negative Binomial Model - Stochastic Block Model
![Negative Binomial Model - Stochastic Block Model](https://github.com/NoemiRossi/Bayesian-Inference-of-Network-Structure-with-Noisy-Data/blob/main/Images/nb_sbm_m.png)
# Poisson Model - Multi-edge Model
![Poisson Model - Multi-edge Model](https://github.com/NoemiRossi/Bayesian-Inference-of-Network-Structure-with-Noisy-Data/blob/main/Images/p_mul_m.png)

# Code and folder structure
* __R codes:__ The folder contains the code that we use to produce our plots and analysis.
* __Datasets:__ The folder contains the dataset used in the Analysis.
* __Stan:__ The folder contains all the .stan file that describes the diffferent models.

# Run a test case
To run a test is enough to run any file in the R codes.

# Results
* __Adjecency Matrix:__From our models we estimate an adjecency matrix of friendship between the students.
![Adjecency Matrix](https://github.com/NoemiRossi/Bayesian-Inference-of-Network-Structure-with-Noisy-Data/blob/main/Images/adj.png)

* __Network:__
![Adjecency Matrix](https://github.com/NoemiRossi/Bayesian-Inference-of-Network-Structure-with-Noisy-Data/blob/main/Images/net.jpg)

# Authors
* __Michele Bellomo__ - *Politecnico di Milano*
* __Paulina Moskwa__ - *Politecnico di Milano*
* __Noemi Rossi__ - *Politecnico di Milano*
