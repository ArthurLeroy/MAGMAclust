This repository contains the R code of the algorithm MAGMAclust, presented in the paper 'Cluster-Specific Predictions with Multi-Task Gaussian Processes' by Leroy, Latouche, Guedj and Gey (https://jmlr.org/papers/v24/20-1321.html). If you are considering citing this work, you are welcome to use the file `reference.bib`.

The folder 'Simulations' contains the synthetic datasets, used for the experimental study in the paper, as well as the corresponding trained models and tables of results.

The file 'Simulation_study.R' contains the code used to generate synthetic datasets, conduct the experiments, evaluate and display the results. 

The file 'Computing_functions.R' contains the code of many useful functions used in MAGMAclust such as kernels, log-likelihoods, gradients, E-step and M-step.

The file 'MAGMAclust.R' contains the core of the code for the MAGMAclust algorithm. Fonctions implementing learning and prediction steps are provided.
Several ploting function and testing examples are provided as well. 

Corresponding author and main contributor to the code:  
Arthur Leroy  
arthur.leroy.pro@gmail.com