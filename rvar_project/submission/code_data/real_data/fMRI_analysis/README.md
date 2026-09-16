# fMRI_analysis

This folder contains the code corresponding to the real-data analysis found in Section 7 in the main paper. First, we here provide a description of the overall structure of the files. Second, we describe in detail the steps required to replicate the simulations, 

## Overall Structure of Files/Directories:

- Directory <code>./req_lib</code>: directory where the required packages will be saved as a local library. It is automatically created when the file <code>./001_requirements.R</code>.

- Files 001-052: these files contain all the required functions to perform simulations. They are divided into specific tasks, such as importing required packages (<code>./001_requirements.R</code>), loading the packages (<code>./002_requirements_lite.R</code>), generating the parameters for performing the real-data analysis on different sample sizes (<code>./051_Simulation_CreatingParameters.R</code>), and evaluate the performance of the fitted models (<code>./052_Simulation_EvaluatingMeasures.R</code>).

- Files 100-110: The files with names 100-110 provide preliminary visualization and processing of the data. In <code>101_CleaningVisualizing_LC.R</code>, we perform preliminary processing, visualize the fMRI time-series and moderator data. In <code>102_AnalysisVarModels.R</code> and <code>103_AnalysisVarModels_1sd.R</code>, we perform model fit on the data for T = 50, 100, and compare the forecasting performance of MOD-VAR and methods in the literature. The RData files  <code>103_Data.RData</code> and <code>110_Data.RData</code> are then used as the processed data for our systematic analysis.

- Files in directories 100-400: These directories contain the calls for the numerical simulations. Each folder corresponds to one method. The files can be used to run (A) a small debugging example to ensure things run smoothly, (B) reduced experiments with a total of 20 simulation replicates, (C) full simulations with 100 simulation replicates.

- Files/Directories 600: After reduced experiments are run, the file <code>./611_Outcomes_DataAggregationRealData.R</code> aggregates the experiments data, and saves it to the directory <code>./600_AggregatedDataReal_v3</code>. Then, file <code>./621_Outcomes_ByPnNSizeExperiment_wFcast.R</code> generates plots with the aggregated data. 

- Files/Directories 700: After full sized numerical experiments, the files <code>711_ModeratorAnalysis.R</code> and <code>712_Plots.R</code> can be used to generate the moderator plots and by-subject decompositions found in Figures 8 and 9 of the main paper.

## Simulation Rerunning Instructions:

The simulations are run in a linux cluster with SLURM scheduling system. 


0. **Clone Repository:** First, clone the github repository <em>MOD-VAR Supplementary Materials</em> to a linux cluster with SLURM scheduling system.

1. **Install Requirements:** Install all packages required for our simulations. For this, in the linux cluster terminal, enter the directory <code>./code_nodata/real_data/fMRI_analysis/</code>. Then run the command

    - <code>Rscript 001_requirements.R</code>

    Ensure that all packages are properly installed by reviewing the directory <code>req_lib/</code>.

4. **Debugging Systematic Simulations:** In the linux cluster terminal, and enter the directory  <code>./simulations/sims_mod2/</code>. To ensure that the code can execute smoothly, run the commands,

    - <code>sbatch 100_var/120_ClusterPassPretrain0.sh</code>
    - <code>sbatch 300_mvar/130_ClusterPassPretrain0.sh</code>
    - <code>sbatch 400_modvar/230_ClusterPassPretrain0.sh</code>
    - <code>sbatch 500_modvar1sd/230_ClusterPassPretrain0.sh</code>
    
    Check the log file <code>100_var/pretrainings1/logs/output0.out</code> to verify that the VAR method ran properly. Often, warnings related to the packages may occur, but check for errors. Check whether the data file <code>100_var/pretrainings1/data/output0_0.RData</code> is created, which ensures that the simulation was completed and successful. You can explore the log and data files in the 200, 300 and 400 directories to ensure the Multi-VAR and MOD-VAR methods also ran successfully. 


3. **Running Simulation Experiments:** Once debugging experiments ran succesfully, you can run preliminary simulation experiments. We have a total of 10 simulation scenarios, corresponding to different choices of T = 10,20,30,...,100. We request 10 cluster nodes, and ask each to perform model fit for each of this sample sizes. To run this, in a cluster terminal, enter the directory <code>./code_nodata/real_data/fMRI_analysis/</code>, and run the lines,

    - <code>sbatch 100_var/130_ClusterPassExperiment0.sh</code>
    - <code>sbatch 300_mvar/130_ClusterPassExperiment0.sh</code>
    - <code>sbatch 400_modvar/230_ClusterPassExperiment0.sh</code>
    - <code>sbatch 500_modvar1sd/230_ClusterPassExperiment0.sh</code>
    
    The logs for the MOD-VAR method experiments are saved in <code>400_modvar/experiments2/logs/</code>, indexed from 1 to 10. The simulation data is saved in RData files in <code>400_modvar/experiments2/data/</code> 1-72. Similar log and data directories for the VAR and Multi-VAR can be found in the folders 100, 300 and 500. 

4. **Generating Simulation Experiment Plots:** Once all preliminary simulation experiments are complete and the data is saved, you can generate simulation plots. For this, in the command line terminal enter the directory <code>./simulations/sims_mod2/</code>, and run

    - <code>Rscript 611_Outcomes_DataAggregationRealData.R</code> 
    
    This saves aggregated data in the directory <code>600_AggregatedDataReal_v3/data_all</code>. Plots that verify performance are saved in the directory <code>600_AggregatedDataReal/plots_all/</code> by running
    
    - <code>Rscript 621_Outcomes_ByPnNSizeExperiment_wFcast.R</code> 
    

6. **Generating Manuscript Plots:** We can generate Figure 7 and 8 of the main paper by running the files: <code>711_ModeratorAnalysis.R</code> and <code>712_Plots.R</code>.
    
