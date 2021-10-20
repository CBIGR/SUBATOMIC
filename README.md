# SUBATOMIC:	SUbgraph BAsed mulTi-OMIcs Clustering
The subgraph based multi-omics clustering (SUBATOMIC) pipeline is a module inference
and annotation framework to inetgrate and cluster networks into functional modules.
Interaction networks can be of any interaction type and contain directed or directed
edges such as TF-target networks, miRNA-target networks or protein-protein interaction
networks as long as they have a common set of shared nodes. 
SUBATOMIC  first integrates all networks into one prior network and decomposes it into 
two- and three-node subgraphs using ISMAGS. The resulting subgraphs are further 
categorized according to their type (ALL, COM, COR, COP, CIR, FBU, FB2U) and
clustered into modules using the hyper edge clustering algorithm SCHype. The resulting
modules contain a high density of subgraphs and are further characterized with 
the functional enrichment method GOATOOLS. Relation with regulators and other
moduels were derived in a superview step summarizing interaction shared between modules and 
regulators pointing at them. 

<object data="http://yoursite.com/the.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="http://yoursite.com/the.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>


## How to run SUBATOMIC: A human mini-net example
Calculating modules for the integrated network published in the paper takes a lot
of compute power and space. We therefore sampled down that network to create 
a mini-net for **H. sapiens** that can be used as a test example. 

### Dependencies
Make sure all dependencies are installed and ready to use. The version number of
python packages should be seen as recommendation, since the pipeline might also 
work with slightly different versions.  

- Python >3.6  
  - numpy 1.20.2  
   - pandas 1.2.4  
   - itertools   
   - subprocess  
   - re  
   - shutil  
   - GOATOOLS 1.0.15  
- Snakemake 6.1.2   
- R > 4.0.0  
   - igraph   

### Installation and start
Clone the repository and adapt the config file.  
The pipeline can be started calling:  
```snakemake Integrated_Pipeline```

### Required and optional input files
Shows an overview of the files required for the analysis. Except for go-basic.obo, all file names can be freely chosen but have to be set in the config file. In this list, all files are contained in the example folder. 

1. ```go-basic.obo```  

- Required for functional enrichment with GOA-tools  
- Located in the base folder next to 'Integrated_Pipeline'
- Can be obtained from [GENEONTOLOGY](http://geneontology.org/docs/download-ontology/)
2. ```goa_human.gaf```  

- Required for functional enrichment with GOA-tools  
- Can be obtained from [GENEONTOLOGY_annotations](http://geneontology.org/docs/download-go-annotations/)

3. ```Hsapens_subgraphs.csv```  

- Plane list of all subgraphs without header
- Due to symmetric reasons, some subgraphs can be redundant. From a set of equivalent subgraphs, only one should be kept
- List can be self defined or generated using the utility script generate_motif_list.py that gives a non-redundant list of motifs given the network letters

4. ```Hsapiens_miniNet.csv```  

- File containing all edges from the integrated prior network
- No header, content:
  - First column contains node_1 (in case of directed interactions the regulators)
  - Second column contains node_2 (in case of directed interactions the targets)
  - Third column contains the letter representing the network to which this interaction belongs
- Networks can only contain **either** directed **or** undirected edges, but no mix of both

5. ```Human_All_Genes.csv```

- Optional, no header, plane list of identifier used as background for functional enrichment analysis
- In this example, a large number of protein coding genes in *H. sapiens*

6. ```Human_functional_description.csv```

- File containing a mapping of identifer to gene name as well as a optional short description of the gene function
- Header, content:
  - First column: Primary identifier, which are also used in the interaction file
  - Second column: External gene name
  - Third column: Potential alternative identifier or gene name
  - Fourth column: functional description

7. ```Human_miRNAS.csv```

- Plane file containing the identifier of all miRNAs
- Will be used in the superview

8. ```Human_TFs.csv```

- Plane file containing the identifier of all regulators
- Will be used in the superview

### Configuration file
 
Check the config file ```config_test_hsapiens_miniNet.yaml```.  
For the toy example, settings are already done.  

- Set all path variables:  
  - **interaction**: *string*,  path + filename to interaction file
  - **motif_type**: *string*, path + filename to file containing all subgraphs to be investigated
  - **Functional_Description**: *string*: path + filename to functinoal description file
  - **Transcription_Factors**: *string*,  path + filename to list of TF
  - **GO_ASSOCIATIONS**: *string*, path + filename to GO annotation file
  - **scripts**: *string*, path to script folder
  - **species**: *string*, any name for an analysis run
  - **directed_interactions**, *list*, letters indicating networks containing exclusively  directed_interactions
  - **undirected_interactions**, *list*, letters indicating networks containing exclusively  undirected_interactions
  - **miRNA_interactions**: *list*, letters indicating the network containing miRNA-target interactions
  - **homolog_interactions**: *list*, letters indicating networks containing undirected homologous interactions  
  - **modules**: *list*, choice from "ALL", "COR", "COP", "CIR", "FBU", "FB2U" and "FFL". Specifies the modules included in the analysis
  - **schype_pvalue**: *int*, weight parameter for SCHYPE (not a p-value actually)
  - **goatools_pvalue**: *float*, threshold p-value for gene enrichment analysis
  -  **call_R**: *string*, command how to call R from the command line. In the easiest case, R
  -  **cores**: *int*, number of cores used for the parallelization of superview. Careful: cores have to be submitted again calling snakemake with --cores. The number of cores specified here needs to be <= the number of --cores.
  - **annotation_background**: *string*, choice what should be taken as background for the functional enrichment ananlysis. Three optinos are available:
      - path + name to a list with gene identifier that will serve as background
      - "ModuleType": All genes of a module type, e.g. all genes contained in a CIR module for CIR
      - "ALL": All genes in the input network
  - **superview_background**: *string*, choice what should be taken as background distribution for the suberview z-score calculation
    - "network": all interactions of a network type
    - "module": all interactions of a network type as well as the module e.g. CIR_H
    
### Running an example analysis for a human miniNet
1. Clone the repository 
2. Go to the base folder and unzip the go-basic.obo file with  
```tar -xzf go-basic.obo.tar.gz``` 
3. Go to example_data/Hsapiens_M1/ and unzip the GAF file with
```tar -xzf goa_human.gaf.tar.gz```
4. Adapt the config file. For this example, all parameters are set and the example is run on 3 cores
5. Go to the SUBATOMIC folder and run
```snakemake -s Integrated_Pipeline --cores 3```
6. Investigate the results

 - All output results can be found in the folder specified in the config file
 - A list and short description of the most important output files can be found
 in the following section 

7. Visualize them with Cytoscpae

- Start Cytoscape
- Open the folder ```/data/human_network/final_networks/Hsapiens/MotifClusters/SCHYPE/SCHype*/```
- The network is contained in the ```edges.nnf file```. Open this file  in Cytoscape (e.g. drag and drop it to the network file area)
- The annotation for this network is included in ```NodeAttributes.noa```.  Open this file  in Cytoscape (e.g. drag and drop it to the table files section)
- The Cytoscape style is included in the folder ```Cytoscape``` and called ```cytoscape_style.xml```. Import the style sheet by opening the ```File``` tab, ```Import``` and ```Styles from File```. This has to be done only once for each combination of directed and un-directed network letters.
- Select all networks, right-click on ```Apply Style``` and select the imported style.  Under the tab Style and Labels, it is possible to change the IDs to gene names. 
- Hint: Cytoscape sometimes does an error in importing naming all listed modules two indices higher (e.g. COR_2 instead of COR_0). In the node table view however, the correct naming is shown. It is also recommended to only load smaller selections of modules at a time to make. 

### Important output files
All results are stored in the designated output folder specified in the description. 
In the following, we describe the content of all folders. However, most files are not
relevant for analysis and interpretation. We therefore first provide a list of most 
important result files:

- ```Cytoscape/cytoscape_style.xml```
- ```GOATOOLS/Enrichment_Modules.csv```
- ```MotifClusters/SCHYPE/SCHype{module_name}/cluster.cla```
- ```MotifClusters/SCHYPE/SCHype{module_name}/edges.nnf```
- ```MotifClusters/SCHYPE/SCHype{module_name}/NodeAttributes.noa```
- ```MotifClusters/Superview/{module type}.csv```
- ```MotifClusters/Superview/miRNA.csv```
- ```MotifClusters/Superview/RF.csv```
- ```MotifClusters/Superview/Module_RF_stats.csv```
- ```MotifClusters/Superview/TF_miRNA_target_stats.csv```

### Overview of all generated output files
A more detailed list of all files can be found here:  
  
1. ```2nodemotifs```  

- Folder
- List of two node motifs found
- Columns:
  - First node
  - Second node
  - Type
- calculated for DD2 (two self-pointing directed interactions)
- calculated for DU2 (two self-pointing interactions, one directed, one undirected)

2. ```3nodemotifs```  

- Folder
- Contains the preprocessed output of ISMAGS for each Motif that was serached for
- First line indicates number of motifs found

3. ```Cytoscape```

- Folder
- Contains the 'cytoscape_style.xml', which is an importable style in cytoscape
- File is created based on the set of undirected and directed nework letters
- Edge colors per network letter can also be changend in this file

4. ```GOATOOLS```

- Folder
- Contains a subfolder per module
- each subfolder contains intermediate results from the parallelized functinoal enrichment analysis
- Important files 
  - ```Enrichment_Modules.csv```
  - Contains all significant functional enticements per module
  - Columns
    - **group**: module type
    - **module_no**: module number
    - **p-value-rank**: per module, what is the term with lowest p-value, second lowest p-value ...
    - **log2foldRatio**: log2 of fold change to indicate strength of enrichment
    - **Remaining**: all remaining columns are default output of GOATOOLS
    
5. ```goatools-main```

- Folder
- Intermediate results for functional enrichment
- Not relevant for final result

6. ```Interactions```

- Folder
- Contains a file for each network type:  ```{species}_{network letter}_{d|u}.txt```
- formatted network files used as input for ISMAGS subgraph search

7. ```MotifClusters```

- Folder
- subfolder ```ModuleViewer``` is not relevant for results 
- subfolder ```SCHYPE``` contains the clustered modules and related information split up for each module type ```SCHype{module_name}```
  - ```cluster.cla```
  - Summary information per cluster
  - Columns
    - **Cluster**: name of the module
    - **NrEdges**: number of edges contained in cluster
    - **PerHOM**: percentage of homologous edges
  -  ```cluster.eda```
  - Assigns each interaction to a cluster
  - Columns
    - First node (e.g. TF)
    - Type of interaction (e.g. R)
    - Second node (e.g. target gene)
    - Assigned cluster
  - ```ClusteringCoefficient.txt```
  - Information about the clustering coefficent per module. A value of 1 indicates a fully connected graph
  - Columns
    - Module name
    - Clustering coefficient
  - ```Data_ALL.txt```
  - Summarizes some more general information per cluster 
  - Columns
    - **Module**: name of the module
    - **Comments**: Annoated comments if some exist
    - **nrTF**: number of TF in the module
    - **nrGO**: number of GO annotated in the module
    - **nrInCluster**: number of genes in the cluster
  - ```edges.nnf```
    - network file for import info Cytoscape
  - ```edges.sif```
    - simple interaction format showing edges and the edge type as extra information
    - Columns
      - source node
      - network type were the edge originates from
      - target node
  - ```NodeAttributes.noa```
  - attribute file that can be imported into Cytoscape to annotate modules
  - information is based on the ```functional_description.csv``` input file
  Columns
    - **Accession_Number**: Gene ID
    - **FuncName**: Name of the gene
    - **Type**: TF, miRNA or gene
    - **TF_type**: family or type of TF
    - **Has_GO**: information whether the gene is annotated with GO-terms
    

8. ```parallel```

- Folder
- Internal folder containing a split of modules into bins for parallelization
- Not relevant as result


9. ```SCHYPE```

- Folder
- Contains the primary output of SCHYpe clustering and some intermediate files
- subfolder ```SCHYPE``` contains the clustered modules and related information split up for each module type ```SCHype{module_name}```
  - ```ALLEDGES.txt```
  -  contains all egdes and their assocaition to a cluster
    - Columns
      - Number of module
      - Source node
      - Edge type, letter indicates network source, upper and lower case indicate edge directions
      - Target node
  - ```*.edges.txt```
  - contain all three node motifs
  - Columns
    - Module number
    - First node
    - Second node
    - Third node
  - ```*.nodes.txt```
    - list for each node to which cluster it belongs
    - Columns
      - Node
      - Module number
  - ```Motifs.txt```
    - list of subgraphs included in all moduels of a certain type
    - Columns
      - First node
      - Second node
      - Third node
  - ```MotifsType.txt```
    - List of motifs and their motif type
    - Columns
      - First node
      - Second node
      - Third node
      - Motif type
      
10. ```Superview```

- Folder
- Contains the result of the superview analysis
- ```{module type}.csv```
  - contain all pais of modules and their interactions
  - Columns
    - **module1**: module type of first module in comparison
    - **module1_nr**: module number of first module in comparison
    - **module1_size**: module size of first module in comparison
    - **module2**: module type of second module in comparison
    - **module2_nr**: module number of second module in comparison
    - **module2_size**: module size of second module in comparison  
    - **{Network_letter}_count**: number of interactions between modules originating from a certain network
    - **{Network_letter}_z-score**: z-score against a comparison of 1000 random modules of the same size 
    - **{Network_letter}_pvalue**: 1-CDF of z-score
  - ```miRNA.csv and RF.csv```
  - summary file of interactions between each regulator and target module
  - Columns
    - **Group**: either miRNA or RF
    - **regulator_name**: identifier of regulator
    - **module2**: module type of  module in comparison
    - **module2_nr**: module number of  module in comparison
    - **module1_size**: module size of  module in comparison
    - **{Network_letter}_count**: number of interactions between RF and module originating from a certain network
    - **{Network_letter}_fraction**: fraction of nodes in the module targeted by a regulator of from a certain network
    - **SumCount**: all counts over all network types
    - **fraction**: fraction of nodes in the module targeted by a regulator     
  - ```Module_RF_stats.csv```
  - Overview per module to summarize how many regulators target this module
  - Columns
    - **Module**: module name
    - **count_RF**: count how many TF target a module
    - **total_RF**: count how many TF were included in the prior network
    - **percentage_RF**: percentage how many RF are targeting a module
    - **count_miRNA**: count how many miRNA target a module
    - **total_miRNA**: count how many miRNA were included in the prior network
    - **percentage_miRNA**: percentage how many RF are targeting a module
  - ```TF_miRNA_target_stats.csv```
  - Overview per module to summarize how many modules are targeted by a specific regulator
  - Columns
    - **RF**: ID of the regulator
    - **type**: type of regulator, either TF or miRNA
    - **count**: number of targeted modules
    - **total**: total number of modules
    - **percentage**: percentage of modules covered comared to the total number of modules      



      
      

 
      



 



