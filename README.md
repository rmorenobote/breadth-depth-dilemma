# breadth-depth-dilemma
Breadth-depth dilemma

This repository contains the cc and Matlab code to carry out the data analysis for the manuscript entitled "Heuristics and optimal solutions to the breadth-depth dilemma" by Rubén Moreno-Bote, Jorge Ramírez-Ruiz, Jan Drugowitsch and Benjamin Y. Hayden.

The cc_code_figures/ folder contains the cc code that produces the data that the Matlab code uses to create each figure. This code was the one used to produce the figures in the manuscript.

The julia_code/ folder contains a Julia version of the code that reproduces the figures in the paper and was used to verify the results. The jupyter notebook "Reproducing_results.ipynb" produces every figure by calling the functions defined in the module "Functions_BD.jl". Here you can find variations of the original code and methods that find the same results.

To run the Julia code, there are two options:
1. If you have Julia v1.0 or higher installed, you can simply download the module "Functions_BD.jl" and the jupyter notebook and run the notebook in the same directory.

2. If you do not have Julia installed, it is very easy to do so (https://julialang.org/downloads/), and install the jupyter notebook (the same that comes with Conda) into it. Another option is to use juliabox (https://juliabox.com/), where you can upload the module and the notebook and run it in the cloud.

