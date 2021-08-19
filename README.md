# Adaptive scoping: balancing short- and long-term genetic gain in plant breeding
## Running the code
- Install the package hypred
- Install the package GSSimTPUpdate
- Open the R project AdaptiveScoping.Rproj
- Run Create_directory.R
- Run MakeGenome_File.R
- Run one of the 'run_experiment' files in the main directory
- Run the correct 'make_figure' script

## make_figures
This directory contains the R script used to make the figures.  

## own_results
Empty directory where simulation results will be saved. This directory is creating by running the script Create_directory.R.

## MakeGenome_File 
Makes a list of n different genomes that can be used in the run_experiment files. At each iteration, each method will use the same genome, making it possible to compare the different methods.

## Genome
Directory containing example files of genomes that can be used in the run_experiment files.

## run_experiment files
R scripts used to simulate a population following the adaptive scoping method or scoping method. The results will be saved in the directory 'own_results'. 
The user will have to load the correct genome from the genome directory and choose the number of nodes that are available for calculation. 

## Supplementary data
Directory containing the supplementary figures and tables.

## R
Directory containing functions that are used in the run_experiment files.