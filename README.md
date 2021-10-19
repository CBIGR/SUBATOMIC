# Composite subgraph cluster integration  
Once we find a better name, we can easily rename the repository.  
Repository is still under construction.  
Starting point is the final version of scripts and pipeline of the Design Project 2020.  

## Content
example_data: contains examples of **D. melanogaster** and **A. thaliana** from the origial approach  
Scripts: contains all binaries and scripts that are called from the snakemake command  
not_integrated_scripts: scripts under construction or scripts not integrated so far in the main pipe  

## How to run
TBA

## Version and hotfix log
All important commits will come with the version number in the commit message.
I will note down in this section what the main changes were.

### Version 0.0.9
- move Logs folder into the dedicated run folder
- delete some commented lines
- repair bug that prevented the superview collector in case there is a module type with no modules in it

### Version 0.0.8
- add pbs scripts to run on HPC
- small modifications and error corrections

### Version 0.0.7
- add utility_scripts folder
- add utility script generate_motif_list.py that allows to generate a non-redundant motif list as ISMAGS input
- add utility script calculate_ECD_score.py, that allows to calculate nPCC and ECD scores. The second is not fully implemented yet, the sampling is not optimized yet
- add functions to new_superview.py, that calculate no. connections between TF and motifs and vice versa to assess specificity

### Version 0.0.6
- Addapting enrichment for parallelsiation
- Add rule to collect parallel enrichment results
- Add script to write a cytoscape style xml that can be used to visualize the results 
- Deprecated Goatools_Significant_enrichment_top5.py

### Version 0.0.5
- Big structural changes introducing a new parallelization scheme
- Adding new rules "parallel", "Super_View_Collector", "Super_View_RF"
- removing rule "improve_SV"
- adding config entries to set the R call and the number of cores one want to use for parallelization
- Adding the parallel_binning.py script to scale up parallelization depending on available nodes
- Re-writing the new_superview.py that merges the superview script functionality with the improve superview R script
- Add new_superview_collector.py that merges results from all parallel runs after they completed
- Add new_superview_RF.py that creates superview files for regulatory factors and miRNAs


### Version 0.0.4
- Adapting SCHYPE to work with user defined network letters
- Removal of calculate_momo in the prepSCHype.py script as well as the consecutive rule
- Smaller adaptions, debugging and improvements to prepSCHype.py (e.g. remoing unnecessary + from w+ arguments in open functions)


### Version 0.0.3
- Integration of ISMAGS into the pipeline
- Adding ISMAGS pre-compiled binary
- Adding script that parses ISMAGS to ISMA output and removes network letters
- update 2node motif detection
- integrate configuration file depended network type integration for all scripts before prepSCHYPE
- removal of ISMA binary and runISMA.py (Move to archive for now)

### Version 0.0.2
- Add a new config file system. Now all configurations and settings will be set in one central file
- Add missing scripts and folders such as post-processing scripts
- Add folder with student feedback (remove in final version)
- Creation of a big bug table showing input, output, descriptions, bugs to fix (not on GitHub, but in Jens Loers' OneDrive)
- Adding archive folder to store scripts not necessary yet but which give valuable insights along the process (good to keep them close)
    
### Version 0.0.1
- Contains the initial version as provided by Design Project 2020


