# simulations


The current directory contains the directories for replicating the systematic numerical simulations presented in Section 6 of our paper, as well as the additional results presented in our Supplementary Materials. In the current README, we provide descriptions for the directories included here. Further descriptions for reproducing our simulation results can be found inside of each of the subdirectories.

- <code>041_modvar/</code>: contains the functions that implement our optimization and tuning procedures for estimating MOD-VAR. These are called by all other subdirectories.

- <code>sims_mod2</code>: contains the code for reproducing our numerical experiments in the case that multi-subject data has heterogeneous GCN, and the heterogeneity is fully explained by the moderators provided. 

- <code>sims_multi2/</code>: contains the code for reproducing our numerical experiments in the case that multi-subject data has heterogeneous GCN, and the heterogeneity is fully individual without moderator effects. This setting is similar to that considered by the Multi-VAR model.

- <code>sims_mpm2/</code>: contains the code for reproducing our numerical experiments in the case that multi-subject data has heterogeneous GCN, and the heterogeneity a combination of moderator effects and fully individual portion.


REMINDER 1: reproduction of the simulation results with the provided code requires clone the directory and all its contents into a linux server with SLURM scheduling system. 

REMINDER 2: Further description of the files and reproducibility can be found in the README.md files found in the directories <code>sims_mod2</code>, <code>sims_multi2/</code> and <code>sims_mpm2/</code>. 