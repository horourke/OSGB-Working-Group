# sims_multi2

This folder contains the code corresponding to the numerical simulations results found in our Supplementary Materials. First, we here provide a description of the overall structure of the files. Second, we describe in detail the steps required to replicate the simulations, 

## Overall Structure of Files/Directories:

- Directory <code>./req_lib</code>: directory where the required packages will be saved as a local library. It is automatically created when the file <code>./001_requirements.R</code>.

- Files 001-003: these files contain all the required functions to perform simulations. They are divided into specific tasks, such as importing required packages (<code>./001_requirements.R</code>), loading the packages (<code>./002_requirements_lite.R</code>), generating synthetic multi-subject time series data (<code>003_Generating_ModVarData.R</code>), generating the parameters for our simulations (<code>./051_Simulation_CreatingParameters.R</code>) and evaluate the performance of the fitted models (<code>./052_Simulation_EvaluatingMeasures.R</code>).

- Files in directories 100-400: These directories contain the calls for the numerical simulations. Each folder corresponds to one method. The files can be used to run (A) a small debugging example to ensure things run smoothly, (B) reduced experiments with a total of 20 simulation replicates, (C) full simulations with 100 simulation replicates.

- Files/Directories 500: After reduced experiments are run, the file <code>./511_Outcomes_DataAggregationExperiments.R</code> aggregates the experiments data, and saves it to the directory <code>./500_AggregatedDataExperiments</code>. Then, file <code>./542_Outcomes_ByPnNSizeExperiment_PAPER.R</code> generates plots with the aggregated data. 

- Files/Directories 600: After full sized simulations with 100 replicates are run, the file <code>./611_Outcomes_DataAggregationFull.R</code> aggregates the full simulation data, and saves it to the folder <code>600_AggregatedDataFull</code>. Then, file <code>642_Outcomes_ByPnNSizeFull_PAPER.R</code> generates plots with the aggregated data. 

## Simulation Rerunning Instructions:

The simulations are run in a linux cluster with SLURM scheduling system. 


0. **Clone Repository:** First, clone the github repository <em>MOD-VAR Supplementary Materials</em> to a linux cluster with SLURM scheduling system.

1. **Install Requirements:** Install all packages required for our simulations. For this, in the linux cluster terminal, enter the directory <code>./simulations/sims_mod2/</code>. Then run the command

    - <code>Rscript 001_requirements.R</code>

    Ensure that all packages are properly installed by reviewing the directory <code>req_lib/</code>.

4. **Debugging Systematic Simulations:** In the linux cluster terminal, and enter the directory  <code>./simulations/sims_mod2/</code>. To ensure that the code can execute smoothly, run the commands,

    - <code>sbatch 100_var/130_ClusterPassExperiment0.sh</code>
    - <code>sbatch 300_mvar/130_ClusterPassExperiment0.sh</code>
    - <code>sbatch 400_modvar/230_ClusterPassExperiment0.sh</code>
    
    Check the log file <code>100_var/experiments1/logs/output0.out</code> to verify that the VAR method ran properly. Often, warnings related to the packages may occur, but check for errors. Check whether the data file <code>100_var/experiments1/data/output0_0.RData</code> is created, which ensures that the simulation was completed and successful. You can explore the log and data files in the 200, 300 and 400 directories to ensure the Multi-VAR and MOD-VAR methods also ran successfully. 


3. **Running Simulation Experiments:** Once debugging experiments ran succesfully, you can run preliminary simulation experiments. We have a total of 72 simulation scenarios, corresponding to different choices of dimension d, number of subjects n, observed time-points per subject T, among other parameters. For each of these simulation scenarios, we perform 20 simulation replicates. To do this, we request 10 cluster nodes, and ask each to perform 2 simulation replicates per simulation scenarios. To run this, in a cluster terminal, enter the directory <code>./simulations/sims_mod2</code>, and run the lines,

    - <code>sbatch 100_var/131_ClusterPassExperiment10.sh</code>
    - <code>sbatch 100_var/132_ClusterPassExperiment20.sh</code>
    - <code>sbatch 300_mvar/131_ClusterPassExperiment10.sh</code>
    - <code>sbatch 300_mvar/132_ClusterPassExperiment20.sh</code>
    - <code>sbatch 400_modvar/231_ClusterPassExperiment10.sh</code>
    - <code>sbatch 400_modvar/232_ClusterPassExperiment20.sh</code>
    
    The logs for the MOD-VAR method experiments are saved in <code>400_modvar/experiments2/logs/</code>, indexed from 10 to 729. The simulation data is saved in RData files in <code>400_modvar/experiments2/data/</code> 1-72. Similar log and data directories for the VAR and Multi-VAR can be found in the folders 100 and 300. 

4. **Generating Simulation Experiment Plots:** Once all preliminary simulation experiments are complete and the data is saved, you can generate simulation plots. For this, in the command line terminal enter the directory <code>./simulations/sims_mod2/</code>, and run

    - <code>Rscript 511_Outcomes_DataAggregationExperiments.R</code> 
    
    This saves aggregated data in the directory <code>500_AggregatedDataExperiments/data_all</code>. Plots that verify performance are saved in the directory <code>500_AggregatedDataExperiments/plots_all/</code> by running
    
    - <code>Rscript 542_Outcomes_ByPnNSizeExperiment_PAPER.R</code> 
    

5. **Running Full Simulations:** You can also run full numerical simulations. For each of the 72 simulation scenarios, we perform 100 simulation replicates. To do this, we request 10 cluster nodes, and ask each to perform 10 simulation replicates per simulation scenarios. To run this, in a cluster terminal, enter the directory <code>./simulations/sims_mod2/</code>, and run the lines,
    
    - <code>sbatch 100_var/141_ClusterPassFull10.sh</code>
    - <code>sbatch 100_var/142_ClusterPassFull20.sh</code>
    - <code>sbatch 300_mvar/141_ClusterPassFull10.sh</code>
    - <code>sbatch 300_mvar/142_ClusterPassFull20.sh</code>
    - <code>sbatch 400_modvar/241_ClusterPassFull10.sh</code>
    - <code>sbatch 400_modvar/242_ClusterPassFull20.sh</code>
    
    The logs for the MOD-VAR full simulations are saved in <code>400_modvar/outputs2/logs/</code>, indexed from 10 to 729. The simulation data is saved in RData files in <code>400_modvar/outputs2/data/</code>. Similar log and data directories for the VAR and Multi-VAR methods can be found in the folders 100 and  300. 

6. **Generating Simulation Experiment Plots:** Once all full simulations are complete, you can generate simulation plots. For this, in the command line terminal enter the directory <code>./simulations/sims_mod2/</code>, and run
    
    - <code>Rscript 611_Outcomes_DataAggregationFull.R</code>
    
    This saves aggregated data in the directory <code>600_AggregatedDataFull/data_all</code>. Plots that verify performance are generated by running
    
    - <code>Rscript 642_Outcomes_ByPnNSizeFull_PAPER.R</code>
    
    Plots are saved in the directory <code>600_AggregatedDataFull/plots_all</code>.