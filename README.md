# breadth-depth-dilemma
Breadth-depth dilemma

This repository contains the Matlab code to carry out the data analysis for the manuscript entitled "Sampling few options is the optimal solution to the breadth-depth dilemma in a finite sample capacity model of decision-making" by Rubén Moreno-Bote, Jorge Ramírez-Ruiz, Jan Drugowitsch and Benjamin Y. Hayden.

The /figures folder contains Matlab code for each of the figures.

The /code folder contains cc code to generate data used in the /figures folder.

The /julia folder contains a Julia version of the code that reproduces the figures in the paper. The jupiter notebook "Reproducing_results.ipynb" produces every figure by calling the functions defined in the module "Functions.jl". This implementation uses the utility estimate that stems from the Markov Chain Monte Carlo method described in the Methods. 

To run the Julia code, there are two options:
1. If you have Julia v1.0 or higher installed, you can simply download the module "Functions.jl" and the jupyter notebook and run the notebook in the same directory.

2. If you do not have Julia installed, it is very easy to do so (https://julialang.org/downloads/), and install the jupyter notebook (the same that comes with Conda) into it. Another option is to use juliabox (https://juliabox.com/), where you can upload the module and the notebook and run it in the cloud.

