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
    parser.add_argument('-run', help='Current run of parallelization')
    parser.add_argument('-groups', help='Modules bzw groups')
    parser.add_argument('-background', help='Choose the background. Options are network or module.')
    parser.add_argument('-interactions', help='Original interaction file.')

    return parser.parse_args()

def intersect(a, b):
    return list(set(a) & set(b))


def readParallelSplit(path):
    """
    Reads in the file containing the split to divide work to the available amount of cores.
    """
    modules = {}
    infile  = open(path, 'r')
    lines   = infile.read().split('\n') 
    
    for line in lines:
        tmp = line.split('_')
        print(tmp)
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
    #if not os.path.exists(f'{Folder}/Superview'):
    #    os.makedirs(f'{Folder}/Superview')

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
    #print(SCHypeGroup + " deleted by sizefilter >= "+str(minSize)+" <= "+str(maxSize)+": " + str(len(deletarray)))
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
    #print(SCHypeGroup + " deleted by homolog filter: " + str(len(deletarray)))
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
    else:
        #print(SCHypeGroup + " detected but not merged clusters with >50% overlap: " + str(len(deletarray)))
        i = 0

    GenesToClusters = {}  # per gene in which cluster they occur
    for Cl in NewNodes.keys():
        Genes = NewNodes[Cl].split('|')
        for G in Genes:
            if G in GenesToClusters.keys():
                GenesToClusters[G] = GenesToClusters[G] + '|' + Cl
            else:
                GenesToClusters[G] = Cl
    
    return NewNodes, GenesToClusters

def superview(int_types, int_data):
    int_out = pd.DataFrame()
    for n in range(0, len(int_data) - 1):
        if (int_data.A[n] in GenesToClusters.keys()) and (int_data.B[n] in GenesToClusters.keys()):  # check if interaction is clustered
            clusters1 = GenesToClusters[int_data.A[n]].split('|')
            clusters2 = GenesToClusters[int_data.B[n]].split('|')
            if clusters1 != clusters2:
                int_out = int_out.append(pd.DataFrame(list(itertools.product(clusters1, clusters2))))
    if int_out.empty:
        return int_out
    int_out.values.sort()  # filter so the smallest clusternumber is fist colum
    int_out = int_out[int_out.apply(lambda x: min(x) != max(x), 1)]  # filter out within modules
    int_count_matrix = int_out.groupby([0, 1]).size().reset_index(name='count')

    int_count_matrix[0] = str(SCHypeGroup + '_') + int_count_matrix[0].astype(str)
    int_count_matrix[1] = str(SCHypeGroup + '_') + int_count_matrix[1].astype(str)
    int_count_matrix['C'] = int_types[0]
    int_count_matrix = int_count_matrix[[0, 'C', 1, 'count']]
    return int_count_matrix

def makePairs(c1, c2):
    interactions = []

    for i in c1:
        for j in c2:
            interactions.append(i+'\t'+j)
            interactions.append(j+'\t'+i)

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

def createSamplesToDrawFromAllInteractions(interactions):
    """
    creates a list of nodes belonging to each random pair (interaction letter, Module)
    e.g 
    """

    new_sets = {}

    infile = open(interactions, 'r')
    lines  = infile.read().split('\n')

    for line in lines:
        tmp = line.split('\t')
        if len(line)>1:
            if tmp[2] not in new_sets:
                new_sets[tmp[2]] = {tmp[0]:'', tmp[1]:''}
            else:
                new_sets[tmp[2]][tmp[0]] = ''
                new_sets[tmp[2]][tmp[1]] = ''

    return new_sets  



def makeSuperview(group, interaction_dict, all_new_nodes, all_gene_to_clusters, lookup_nodes, split, path, groups, background, interactions):
    """
    Create the superview (interactions between clusters)
    """    
    
    # get a sampling pool
    if background == "module":   
        pool = createSamplesToDrawFrom(lookup_nodes, all_new_nodes, groups)     
    else:
        pool = createSamplesToDrawFromAllInteractions(interactions)
        
    # get network file letters
    network_files = interaction_dict.keys()

    # storage of random distributions
    storage = {}

    # prepare file to write
    
    if os.path.exists(path):
        outfile = open(path, 'a+')
    else:
        outfile = open(path, 'w+')
        header  = "module1\tmodule1_nr\tmodule1_size\tmodule2\tmodule2_nr\tmodule2_size\t"
        for i in network_files:
            header += i+'_'+'count'+"\t"+i+'_'+'zscore'+"\t"+i+'_'+'pvalue'+"\t"
        outfile.write(header[:-1]+'\n')    

    # go through all modules and calculate their interaction with other modules

    for i in all_new_nodes[group]:

        # skip if i is not in the parallel split
        if i not in split:
            continue
        
        cluster1 = all_new_nodes[group][i].split('|')

        # compare with all other clusters
        for group2 in groups:
            for k in all_new_nodes[group2]:
                cluster2 = all_new_nodes[group2][k].split('|')
                
                # make all pairs between those clusters                
                pairs = makePairs(cluster1, cluster2)
                # initialize count dict storing interactions between moduels
                count = {}
                for c in network_files:
                    count[c] = 0
                count['all'] = 0

                # check existence of pairs in the interaction file
                for pair in pairs:
                    for n in network_files:
                        if pair in interaction_dict[n]:
                            count[n] += 1
                            count['all'] += 1
                #print(group,group2, i,k, len(cluster1), len(cluster2))
                #print(count)

                # From here, one can calculate the improved superview
                string = [group, i,len(cluster1), group2, k,len(cluster2)]
                #print(count)
                for n in network_files:
                    #string.append(n)
                    marker = n+'_'+group+'_'+str(len(cluster1))+'_'+group2+'_'+str(len(cluster2))

                    if marker in storage:
                        if count[n] != 0:
                            z_score  = (count[n] - storage[marker][0] ) / (storage[marker][1]+0.001)
                            p_values = 1- stats.norm.cdf(abs(z_score))
                            if z_score < 0:
                                z_score  = 0
                                p_values = 1
                        else:
                            z_score  = 0
                            p_values = 1
                        string.append(count[n])
                        string.append(round(z_score,2))
                        string.append(p_values)
                    else:
                        distribution = []
                        
                        if count[n] != 0:
                            warning = False
                            for c in range(1,1000):
                                # in case of module sample for network and module type (n+group). In other case, sample from the complete space of n
                                if background == "module":
                                    if n+'_'+group in pool:
                                        # extra condition to check size of sampling pool. If there is a problem, write a warning
                                        # This can happen if interactions appear twice (e.g. R is also H). Then also for a CIR, the
                                        # H is counted. This cannot be prevented in this current version, since that information is 
                                        # lost in the SCHYPE clustering
                                        if len(cluster1) <= len(pool[n+'_'+group]) and len(cluster2) <= len(pool[n+'_'+group2]):
                                            set1 =  random.sample(pool[n+'_'+group].keys(),len(cluster1))
                                            set2 =  random.sample(pool[n+'_'+group2].keys(),len(cluster2))
                                        else:
                                            #print("Warning: insufficient sampling pool exist for "+ n+'_'+group+i+'_'+group2+k+'. P-value set to 1')
                                            set1 = []
                                            set2 = []
                                            warning = True
                                    else:
                                        set1 = []
                                        set2 = []
                                else:
                                    if n in pool:
                                        # extra condition to check size of sampling pool. If there is a problem, write a warning
                                        # This can happen if interactions appear twice (e.g. R is also H). Then also for a CIR, the
                                        # H is counted. This cannot be prevented in this current version, since that information is 
                                        # lost in the SCHYPE clustering
                                        
                                        if len(cluster1) <= len(pool[n]) and len(cluster2) <= len(pool[n]):
                                            set1 =  random.sample(pool[n].keys(),len(cluster1))
                                            set2 =  random.sample(pool[n].keys(),len(cluster2))
                                        else:
                                            #print("Warning: insufficient sampling pool exist for "+ n+'_'+group+i+'_'+group2+k+'. P-value set to 1')
                                          #print(len(cluster1), len(pool[n]), len(cluster2), len(pool[n]))
                                            set1 = []
                                            set2 = []
                                            warning = True
                                    else:
                                        set1 = []
                                        set2 = []
                                #print(distribution)         
                                rand_pairs = makePairs(set1, set2)
                                interaction_counter = 0
                                for pair in rand_pairs:
                                    if pair in interaction_dict[n]:
                                        interaction_counter +=1
                                distribution.append(interaction_counter)
                            # print a warning in case a certain cluster has not enough info
                            if warning != True:
                                #
                                
                                                          
                                # store values for later use against same distributions
                                storage[marker] = [round(np.mean(distribution),3), round(np.std(distribution),3)]
                                # calculate z-score
                                z_score  = (count[n] - round(np.mean(distribution),3) ) / (round(np.std(distribution),3)+0.001)
                                p_values = 1- stats.norm.cdf(abs(z_score))

                                # remove negative values
                                if z_score < 0:
                                    z_score  = 0
                                    p_values = 1

                                string.append(count[n])
                                string.append(round(z_score,2))
                                string.append(p_values)

                            else:
                                print("Warning: insufficient sampling pool exist for "+ n+'_'+group+i+'_'+group2+k+'. P-value set to NA')
                                z_score  = 'NA'
                                p_values = 'NA'
                                string.append(count[n])
                                string.append(z_score)
                                string.append(p_values)
                                warning = False
                        else:
                            z_score  = 0
                            p_values = 1
                            string.append(count[n])
                            string.append(round(z_score, 2))
                            string.append(p_values)
                


                outfile.write("\t".join([str(item) for item in string])+'\n')
                    
       

if __name__ == '__main__':
    """
    Script that combines the creation of a superview with its statistical evaluation.
    It is supposed to work per module type and needs to be as computational efficient as possible
    It tries to save as much code (and work) from orignial proposal
    Possible to select between the TF mode and the Regulator mode
    """

    #parse argumnets
    args = arg_parse()
    

    # load split for parallelization
    split = readParallelSplit(args.run)

    # groups to run against
    groups = args.groups.split('_')
    
    # calculate interactions between all modules
    if args.mode == 'module':
        # prepare folder and read in all the stuff as dict
        interaction_dict, lookup_ints, lookup_nodes = loadInteractions(args.folder)
        #exit()
        all_new_nodes = {}
        all_gene_to_clusters = {}
        # load modules
        for SCHypeGroup in groups:#Groups:
            NewNodes, GenesToClusters = loadModules(args.folder, SCHypeGroup, 5, 50, args.hom, False)
            print(SCHypeGroup, ": ", len(NewNodes), " found", )
            all_new_nodes[SCHypeGroup] = NewNodes
            all_gene_to_clusters[SCHypeGroup] = GenesToClusters

        
        ## calculate superview
        name = args.run.replace("\\","/").split('/')[-1]
        path = args.folder+'/Superview/'+name
        if os.path.exists(path):
            os.remove(path)

        for SCHypeGroup in split.keys():
            makeSuperview(SCHypeGroup, lookup_ints, all_new_nodes, all_gene_to_clusters, lookup_nodes, split[SCHypeGroup], path, groups, args.background, args.interactions)

    # calculate interactions between TF and all modules
    if args.mode == 'rf':
        exit()

#Logs/Super_View_{species}_{runs}_done.txt
    name = args.run.replace("\\","/").split('/')[-1]

    try:
        os.mkdir(args.folder+"/Logs")
    except FileExistsError:
        pass
    with open(args.folder+"/Logs/Super_View_"+args.folder+"_"+name.split('.')[0]+"_done.txt", 'w+') as f:
        f.write("Super view done!")




    
