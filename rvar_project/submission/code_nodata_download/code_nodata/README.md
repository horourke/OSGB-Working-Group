# Mod-VarCodeSubmission
This repository contains all the proof-read and commented code to recreate simulations, figures and tables of the paper "Quantifying Moderator Effects on High-Dimensional Multi-Subject Time Series"

Here, we simply describe the overall structure of the repository. To fully reproduce the results, further instructions are necessary. For this end, we provide further README files when necessary. The structure is the following:

- <code>./figures/</code>: directory containing the files with all figures generated manually, and do not require extensive simulations or real data. 

- <code>./real_data/</code>: directory containing the code for performing the real data analysis we describe in Section 7 of our main paper. Due to privacy concerns, we do not provide the raw real data, but we provide the processed results, and a README file providing further instructions on how to replicate the production of the figures from the real data analysis section. 

- <code>./simulations/</code>: directory containing the code for reproducing the simulation results of Section 6 in our main paper, and Section XXX of our Supplementary Materials. In order to replicate all simulation experiments and plots, you must run the simulation code on a linux cluster with SLURM scheduling system. We provide further description of how to reproduce the results in the README of this directory.