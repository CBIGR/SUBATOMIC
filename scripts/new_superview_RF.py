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
    parser.add_argument('-mode', help='Mode of integration. Can be Regualtory factor (rf) or (module) .')
    parser.add_argument('-groups', help='Modules bzw groups')
    parser.add_argument('-RF', help='Path to file holding regulatory factors')
    parser.add_argument('-miRNA', help='Path to file holding miRNAs')
    #parser.add_argument('Path_to_SCHYPE_jar', help='The relative or absolute path to the master-java folder containing the SCHYPE java script.')

    return parser.parse_args()

def intersect(a, b):
    return list(set(a) & set(b))

def loadRF(path):
    """
    Load the TF file
    """
    TF = {}
    infile  = open(path, 'r')
    lines   = infile.read().split('\n') 
    
    for line in lines:
        tmp = line.split('\t')
        if len(line) > 0:
            TF[tmp[0]] = ''
    
    return(TF)                

def readParallelSplit(path):
    """
    Reads in the file containing the split to divide work to the available amount of cores.
    """
    modules = {}
    infile  = open(path, 'r')
    lines   = infile.read().split('\n') 
    
    for line in lines:
        tmp = line.split('_')
        if len(tmp) > 1:
            if tmp[0] in modules:
                modules[tmp[0]][tmp[1]] = ''
            else:
                modules[tmp[0]] = {tmp[1]:''}
    return(modules)                



def loadInteractions(Folder):
    """
    Ceate all superview folders and afterwards create a dictionary containing all interactions 
    """
    #Groups = [G.replace('SCHype', '') for G in os.listdir(f'{Folder}/SCHYPE')]  # UNDEF, MOMO, CIR, COR, COP, COM, FBU, FB2U, FFL, ALL
    if not os.path.exists(f'{Folder}/Superview'):
        os.makedirs(f'{Folder}/Superview')

    interactions = [f'{Folder}/Interactions/{G}' for G in os.listdir(f'{Folder}/Interactions')]
    interaction_dict = {}
    for data in interactions:
        if data[-7].isupper():
            interaction_dict[data[-7] + 'int'] = pd.read_csv(data, sep='\t', header=0, names=['A', 'B'])
    
    # make interaction dictionary that can be efficiently used for superview
    ints = {}
    nodes = {}
    for data in interactions:
        d1 = {}
        nodes1 = {} 
        infile = open(data, 'r')
        lines = infile.read().split('\n')

        for line in lines:
            if len(line) >0:
                tmp = line.split('\t')
                d1[line] = ''
                nodes1[tmp[0]] = ''
                nodes1[tmp[1]] = ''
        
        ints[data[-7]] = d1
        nodes[data[-7]] = nodes1
    return interaction_dict, ints, nodes

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
        for ElementToDelete in deletarray:
            del NewNodes[ElementToDelete]
    else:
        u=0

    GenesToClusters = {}  # per gene in which cluster they occur
    for Cl in NewNodes.keys():
        Genes = NewNodes[Cl].split('|')
        for G in Genes:
            if G in GenesToClusters.keys():
                GenesToClusters[G] = GenesToClusters[G] + '|' + Cl
            else:
                GenesToClusters[G] = Cl
    
    return NewNodes, GenesToClusters


def makePairs(c1, c2):
    interactions = []

    for i in c1:
        for j in c2:
            interactions.append(i+'\t'+j)

    return interactions

def createSamplesToDrawFrom(lookup_nodes, all_new_nodes, groups):
    """
    creates a list of nodes belonging to each random pair (interaction letter, Module)
    e.g 
    """

    new_sets = {}

    for n in lookup_nodes:

        for group in groups:
            mini_set = {}
            for i in all_new_nodes[group]:
                cluster1 = all_new_nodes[group][i].split('|')
                for i in cluster1:
                    if i in lookup_nodes[n]:
                        mini_set[i] = ''
            new_sets[n+'_'+group] = mini_set
        
    for group in groups:
        mini_set = {}
        for i in all_new_nodes[group]:
            cluster1 = all_new_nodes[group][i].split('|')             
            for i in cluster1:
                mini_set[i] = ''
        new_sets['all_'+group] = mini_set
      
    return new_sets        



def makeSuperview(group, interaction_dict, all_new_nodes, all_gene_to_clusters, lookup_nodes, groups, path, factors):
    """
    Create the superview (interactions between clusters)
    """    
    

    # get network file letters
    network_files = interaction_dict.keys()

    # store information for further 
    all_interactions = []
    # prepare file to write
   
    outfile = open(path, 'w+')
    header  = "Group\tregulator_name\tmodule2\tmodule2_nr\tmodule2_size\t"
    for i in network_files:
        header += i+'_'+'count'+"\t"+i+'_'+'fraction'+"\t"
    header += "SumCount\tfraction"    
    outfile.write(header+'\n')  
    
    if factors == None or factors == []:
        return 

    # go through all modules and calculate their interaction with other modules
    for TF in factors:
        #print(TF)
        # compare with all other clusters
        for group2 in groups:
            for k in all_new_nodes[group2]:
                cluster2 = all_new_nodes[group2][k].split('|')
                # make all pairs between those clusters                
                pairs = makePairs([TF], cluster2)
                # initialize count dict storing interactions between moduels
                count = {}
                for c in network_files:
                    count[c] = 0
                
                all_counts = {}
                # check existence of pairs in the interaction file
                for pair in pairs:
                    for n in network_files:
                        if pair in interaction_dict[n]:
                            count[n] += 1
                            all_counts[pair] = 0
                
                string = [group, TF, group2, k,len(cluster2)]
                for c in count:
                    string.append(count[c])
                    string.append(round(count[c]/len(cluster2),4)*100)
                string.append(len(all_counts))
                string.append(round(len(all_counts)/len(cluster2),4)*100)
                outfile.write("\t".join([str(item) for item in string])+'\n')
                all_interactions.append(string)

    # calcaute How many TF target a specitic 
    return all_interactions
                  

def counter(TF_interactions, miRNA_interactions, path):
    """
    Make statistic of how many TF target a specific module and how many modules are targeted by one TF
    """

    RF_Counter        = {}
    Module_counter_mi = {}
    Module_counter_tf = {}
    types             = {}
    module_count      = {}
    TFs               = {}


    for i in TF_interactions:
        if i[1] in RF_Counter and i[-2]>0:
            RF_Counter[i[1]][i[2]+'_'+i[3]] = ''                
        else:
            if i[-2]>0:
                RF_Counter[i[1]] = {i[2]+'_'+i[3]:''}
                types[i[1]] = i[0]

        if i[2]+'_'+i[3] in Module_counter_tf and i[-2]>0:
            Module_counter_tf[i[2]+'_'+i[3]][i[1]] = ''
        else:
            if i[-2]>0:
                Module_counter_tf[i[2]+'_'+i[3]] = {i[1]:''}
                types[i[1]] = i[0]
                module_count[i[2]+'_'+i[3]] = ''

    for i in miRNA_interactions:
        if i[1] in RF_Counter and i[-2]>0:
            RF_Counter[i[1]][i[2]+'_'+i[3]] = ''                
        else:
            if i[-2]>0:
                RF_Counter[i[1]] = {i[2]+'_'+i[3]:''}
                types[i[1]] = i[0] 

        if i[2]+'_'+i[3] in Module_counter_mi and i[-2]>0:
            Module_counter_mi[i[2]+'_'+i[3]][i[1]] = ''
        else:
            if i[-2]>0:
                Module_counter_mi[i[2]+'_'+i[3]] = {i[1]:''}
                types[i[1]] = i[0] 
                module_count[i[2]+'_'+i[3]] = ''

    outfile1 = open(path+'/TF_miRNA_target_stats.csv', 'w+')
    outfile2 = open(path+'/Module_RF_stats.csv', 'w+')
    header1  = "RF\ttype\tcount\ttotal\tpercentage\n"
    header2  = "Module\tcount_RF\ttotal_RF\tpercentage_RF\tcount_miRNA\ttotal_miRNA\tpercentage_miRNA\n"
    outfile1.write(header1)
    outfile2.write(header2)

    for i in RF_Counter:
        outfile1.write('\t'.join([i, types[i], str(len(RF_Counter[i])),str(len(module_count)), str(round(len(RF_Counter[i])/len(module_count),4)*100), '\n']))

  
    lenRF = 0
    lenMI = 0
    
    for i in types:
        if types[i] == 'TF':
            lenRF += 1
        else:
            lenMI += 1

    bl = {}

    for i in Module_counter_tf:
        
        if i not in Module_counter_mi:
            outfile2.write('\t'.join([i, str(len(Module_counter_tf[i])),str(lenRF), str(round(len(Module_counter_tf[i])/lenRF,4)*100),'0','0','0', '\n']))
        else:
            outfile2.write('\t'.join([i, str(len(Module_counter_tf[i])),str(lenRF), str(round(len(Module_counter_tf[i])/lenRF,4)*100), str(len(Module_counter_mi[i])),str(lenMI),str(round(len(Module_counter_mi[i])/lenRF,4)*100), '\n']))
            bl[i] = ''
                  
    for i in Module_counter_mi:
        
        if i not in bl:
            outfile2.write('\t'.join([i,'0','0','0',str(len(Module_counter_mi[i])),str(lenRF), str(round(len(Module_counter_mi[i])/lenRF,4)*100), '\n']))
            
                  

       

if __name__ == '__main__':
    """
    Script that combines the creation of a superview for regulatory factors TF and miRNA
    It is supposed to work per module type and needs to be as computational efficient as possible
    It tries to save as much code (and work) from orignial proposal
    Possible to select between the TF mode and the Regulator mode
    """

    #parse argumnets
    args = arg_parse()
    
    # groups to run against
    groups = args.groups.split('_')
    
    # load TF
    if os.path.exists(args.RF):
        TF = loadRF(args.RF)
    else:
        TF = None

    # load miRNA
    if os.path.exists(args.miRNA):
        miRNA = loadRF(args.miRNA)
    else:
        miRNA = None
    
    # calculate interactions between all modules
    if args.mode == 'rf':
        # prepare folder and read in all the stuff as dict
        interaction_dict, lookup_ints, lookup_nodes = loadInteractions(args.folder)
        #exit()
        all_new_nodes        = {}
        all_gene_to_clusters = {}
        # load modules
        for SCHypeGroup in groups:#Groups:
            NewNodes, GenesToClusters = loadModules(args.folder, SCHypeGroup, 5, 50, args.hom, False)
            print(SCHypeGroup, ": ", len(NewNodes), " found", )
            all_new_nodes[SCHypeGroup] = NewNodes
            all_gene_to_clusters[SCHypeGroup] = GenesToClusters

        
        ## calculate superview

        # for TF
        path = args.folder+'/Superview/RF.csv'
        TF_interactions = makeSuperview("TF", lookup_ints, all_new_nodes, all_gene_to_clusters, lookup_nodes, groups, path, TF)

        # for miRNA
        path = args.folder+'/Superview/miRNA.csv'
        miRNA_interactions = makeSuperview("miRNA", lookup_ints, all_new_nodes, all_gene_to_clusters, lookup_nodes, groups, path, miRNA)
        

        # calculate How many TF or miRNA target a module and vice versa 
        counter(TF_interactions, miRNA_interactions, args.folder+'/Superview/')

    try:
        os.mkdir(args.folder+"/Logs")
    except FileExistsError:
        pass
    with open(args.folder+"/Logs/Super_View_RF_done.txt", 'w+') as f:
        f.write("Super view done!")




    
