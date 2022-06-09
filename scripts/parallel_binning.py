import sys
import os
import csv
import pandas as pd
import itertools
import os
import sys
import argparse
import random
import numpy as np
from scipy import stats
import time

def arg_parse():
    """
    Parses the arguments given in the command line
    """
    parser = argparse.ArgumentParser(description='Prep and execute SCHYPE.')
    parser.add_argument('-folder', help='The relative or absolute path to folder where your input files are located, in this folder, also the output will be written.')
    parser.add_argument('-hom', help='List of homolog network letters')
    #parser.add_argument('-cores', help='List the undirected interaction letters.')
    #parser.add_argument('Path_to_SCHYPE_jar', help='The relative or absolute path to the master-java folder containing the SCHYPE java script.')

    return parser.parse_args

def intersect(a, b):
    return list(set(a) & set(b))

def loadModules(Folder, Groups, minSize, maxSize, hom, merge):
    """
    Load and filter schype modules. Output:
    NewNodes: contains all modules and its genes
    
    """

    # translate modules to one node with gene list as feature
    NewNodes = {}
    with open((Folder + '/SCHYPE/SCHype' + SCHypeGroup + '/' + SCHypeGroup + '.nodes.txt'), 'r') as nodefile:
        nodescluster = csv.DictReader(nodefile, delimiter='\t', dialect="excel", fieldnames=['gene', 'cluster'])
        for row in nodescluster:
            if row['cluster'] in NewNodes.keys():
                NewNodes[row['cluster']] = NewNodes[row['cluster']] + '|' + row['gene']
            else:
                NewNodes[row['cluster']] = row['gene']
    
    
    # filter cluster nodes between minSize and maxSize
    deletarray = []
    ClusterSize = {}
    for Cl in NewNodes.keys():
        Genes = NewNodes[Cl].split('|')
        ClusterSize[Cl] = len(Genes)
        if len(Genes) < minSize:
            deletarray.append(Cl)
        elif len(Genes) > maxSize:
            deletarray.append(Cl)

    # delete small/big cluster
    for ElementToDelete in deletarray:
        del NewNodes[ElementToDelete]


    # Delete full homolog clusters
    Hcount = {}
    EdgeCount = {}
    with open((Folder + '/SCHYPE/SCHype' + SCHypeGroup + '/ALLEDGES.txt'), 'r') as edgefile:
        edgecluster = csv.DictReader(edgefile, delimiter='\t', dialect="excel",
                                     fieldnames=['cluster', 'gene1', 'interaction', 'gene2'])
        for row in edgecluster:
            if row['cluster'] in Hcount.keys() and row['interaction'] in hom:
                Hcount[row['cluster']] = Hcount[row['cluster']] + 1
            elif row['interaction'] in hom:
                Hcount[row['cluster']] = 1

            if row['cluster'] in EdgeCount.keys():
                EdgeCount[row['cluster']] = EdgeCount[row['cluster']] + 1
            else:
                EdgeCount[row['cluster']] = 1

    deletarray = []
    for Cl in Hcount.keys():
        if Hcount[Cl] / EdgeCount[Cl] > 0.9:
            deletarray.append(Cl)

    # delete homologs
    for ElementToDelete in deletarray:
        if ElementToDelete in NewNodes.keys():
            del NewNodes[ElementToDelete]

    

    # Merge clusters which are 50 % overlapping


    deletarray = []
    for Cl in NewNodes.keys():
        Genes = NewNodes[Cl].split('|')
        for Cl2 in NewNodes.keys():
            if Cl != Cl2:
                Genes2 = NewNodes[Cl2].split('|')
                if len(intersect(Genes, Genes2)) / len(list(set(Genes + Genes2))) > 0.5:

                    # merge small one in bigger one, delete te small one
                    if len(Genes) >= len(Genes2):
                        NewNodes[Cl] = '|'.join(list(set(Genes + Genes2)))
                        if Cl2 not in deletarray:
                            deletarray.append(Cl2)
                    else:
                        NewNodes[Cl2] = '|'.join(list(set(Genes + Genes2)))
                        if Cl not in deletarray:
                            deletarray.append(Cl)

    

    if merge == True:
        #print(SCHypeGroup + " detected and merged clusters with >50% overlap: " + str(len(deletarray)))
        for ElementToDelete in deletarray:
            del NewNodes[ElementToDelete]
    #else:
        #print(SCHypeGroup + " detected but not merged clusters with >50% overlap: " + str(len(deletarray)))


    GenesToClusters = {}  # per gene in which cluster they occur
    for Cl in NewNodes.keys():
        Genes = NewNodes[Cl].split('|')
        for G in Genes:
            if G in GenesToClusters.keys():
                GenesToClusters[G] = GenesToClusters[G] + '|' + Cl
            else:
                GenesToClusters[G] = Cl
    
    return NewNodes, GenesToClusters




if __name__ == '__main__':
    """
    Create bins that are can be used for parallelization
    """

    #args = arg_parse()
    folder = sys.argv[1]
    hom    = sys.argv[2]
    cores  = int(sys.argv[3])
    modules = sys.argv[4].split('_')

   
    if not os.path.exists(f'{folder}/parallel'):
        os.makedirs(f'{folder}/parallel')

    # additionally, create the Superview folder already here:
    if not os.path.exists(f'{folder}/Superview'):
        os.makedirs(f'{folder}/Superview')


#    else:
#        af = os.listdir(f'{folder}/parallel/')
#        for i in af:        
#            os.remove(f'{folder}/parallel/'+i)

    all_new_nodes = {}
    all_gene_to_clusters = {}
    # load modules
    summer = 0
    for SCHypeGroup in modules:#Groups:
        NewNodes, GenesToClusters = loadModules(folder, SCHypeGroup, 5, 50, hom, False)
       # print(SCHypeGroup, ": ", len(NewNodes), " found", )
        summer += len(NewNodes)
        all_new_nodes[SCHypeGroup] = NewNodes
        all_gene_to_clusters[SCHypeGroup] = GenesToClusters

    ticks = int(summer/cores)+1
    counter  = 0
    counter2 = 1

    outfile = open(f'{folder}/parallel/0.txt', 'w')

    for SCHypeGroup1 in modules:
    #for SCHypeGroup1 in ["COM"]:
        
        for cluster1 in all_new_nodes[SCHypeGroup1]:
            #print(counter, SCHypeGroup1, cluster1)
            outfile.write(SCHypeGroup1+'_'+cluster1+'\n')
            if counter2 != ticks:
                counter2 +=1
            else:
                counter += 1
                counter2 = 1
                outfile.close()
                outfile = open(f'{folder}/parallel/'+str(counter)+'.txt', 'w')
                
        #counter +=1 
                

    #print(int(summer/cores)+1)

                                                  
                                

    
    #for i in range(0,)


   


try:
    os.mkdir(f"{folder}/Logs")
except FileExistsError:
    pass
with open(f"{folder}/Logs/parallel_done.txt", 'w') as f:
    f.write("Super view done!")


