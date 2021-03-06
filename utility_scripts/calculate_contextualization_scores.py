import argparse
from os import listdir
from os.path import isfile, join
import os
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
import warnings
warnings.filterwarnings("ignore")

def arg_parse():
    """
    Parses the arguments given in the command line
    """
    parser = argparse.ArgumentParser(description='Calculates ECD scores, nPCC and biological activity per module given a set of expression data or a list of p-values from differential expression.')
    parser.add_argument('-expression', help='Path + name of expression file')
    parser.add_argument('-pvalues', help='Path + name of file containing differentialy expressed pvalues')
    parser.add_argument('-meta', help='Path + name meta data table')
    parser.add_argument('-clusters', help='Path to the SCHYPE output folder (containing SCHypeALL, SCHypeCIR etc.') 
    parser.add_argument('-out', help='Path to the folder to store results')
    parser.add_argument('-nPCC', action="store_true", help='Flag to include nPCC calculation')
    parser.add_argument('-ECD', action="store_true", help='Flag to include ECD calculation')
    parser.add_argument('-activity', action="store_true", help='Flag to include ECD calculation')
    parser.add_argument('-sampling_k', help='Number of samples for randomization')

    return parser.parse_args()

def loadMeta(path):
    """
    Load the meta data file (right now only interesting spliting up the header into case and controll)
    """
    
    infile = open(path, 'r')
    lines  = infile.read().split('\n')

    case    = {}
    control = {}
    experiments = {}


    # load full data into the experiments folder
    for line in lines[1:]:
        tmp = line.split('\t')
        
        if len(tmp) > 1:
            if tmp[3] == '0':
                # create experiment entry if not existing
                if tmp[1] not in experiments:
                   experiments[tmp[1]] = {}
                # add condition if it exisist
                if tmp[2] in experiments[tmp[1]]:
                   experiments[tmp[1]][tmp[2]].append(tmp[0])
                else:
                   experiments[tmp[1]][tmp[2]] = [tmp[0]]
                case[tmp[0]] = ''
                    
            
            if tmp[3] == '1':
                # create experiment entry if not existing
                if tmp[1] not in experiments:
                   experiments[tmp[1]] = {}
                # add condition if it exisist
                if 'control' in experiments[tmp[1]]:
                   experiments[tmp[1]]['control'].append(tmp[0])
                else:
                   experiments[tmp[1]]['control'] = [tmp[0]]
            control[tmp[0]] = ''
    return experiments, case, control

def pooling(modules):
    """
    Merge all module genes into one pool
    """
    
    pool = {}

    for module in modules:
        for gene in modules[module]:
            pool[gene] = ''
    return pool.keys()
    

def loadExpression(path):
    """
    Load expression data from csv and store it in a dictionary per gene
    """

    infile = open(path, 'r')
    lines  = infile.read().split('\n')

    genes  = {}
    header = lines[0].split('\t')[1:]

    for line in lines[1:]:
        tmp = line.split('\t')
        if len(tmp) > 1:
            genes[tmp[0]] = tmp[1:]

    return header, genes

def loadPvalues(path):
    """
    Load p-values from differential expression file per gene
    """

    infile = open(path, 'r')
    lines  = infile.read().split('\n')

    pvalues  = {}

    for line in lines[1:]:
        tmp = line.split('\t')
        if len(tmp) > 1:
            if tmp[0] not in pvalues:
                pvalues[tmp[0]] = float(tmp[1])
            else:
                if pvalues[tmp[0]] < float(tmp[1]):
                    pvalues[tmp[0]] = float(tmp[1])
        

    return pvalues
         
def loadModules(path):
    """
    Load the modules generated by SCHYPE
    """
    
    path_content = listdir(path)

    modules  = {}

    for i in path_content:
        motif_type = i[6:]
        infile     = open(path+'/'+i+'/'+motif_type+'.nodes.txt', 'r')
        lines      = infile.read().split('\n')

        for line in lines:
            tmp = line.split('\t')
            
            if len(tmp) > 1:
                if motif_type+'_'+tmp[1] in modules:
                    modules[motif_type+'_'+tmp[1]].append(tmp[0])
                else:
                    modules[motif_type+'_'+tmp[1]] = [tmp[0]]
    return modules

def loadInteractions(path):
    """
    Load the interactions generated from ALLEDGES. Makes an entry for 
    both directions since we don't know in which order they are requested later
    """
    
    path_content = listdir(path)

    modules  = {}

    for i in path_content:
        motif_type = i[6:]
        infile     = open(path+'/'+i+'/ALLEDGES.txt', 'r')
        lines      = infile.read().split('\n')

        for line in lines:
            tmp = line.split('\t')
            
            if len(tmp) > 1:
                if tmp[1]+'_'+tmp[3] in modules:
                    modules[tmp[1]+'_'+tmp[3]].append(tmp[0])
                else:
                    modules[tmp[1]+'_'+tmp[3]] = [tmp[0]]
                if tmp[3]+'_'+tmp[1] in modules:
                    modules[tmp[3]+'_'+tmp[1]].append(tmp[0])
                else:
                    modules[tmp[3]+'_'+tmp[1]] = [tmp[0]]
    return modules

def make_artifical_module(list_of_items, l):
    new_list = {}
    
    for i in list_of_items:
        tmp = i.split('_')
        if len(new_list) < l:
            new_list[tmp[0]] = ''
        if len(new_list) < l:
            new_list[tmp[1]] = ''
            
    return list(new_list.keys())
    

def ECDscore(outfolder, experiments, case, control, header, genes, modules, interactions, k):
    """
    Calculate ECD score per module.
    """

    # create Folder to store 
    if os.path.isdir(outfolder+'/ECD') == False:
        os.mkdir(outfolder+'/ECD')

    outfile = open(outfolder+'/ECD/ECD.csv', 'w')
    outfile.write("module\texperiment\tcondition\tlength_module\tmean_dPCC\tz_score\tp_values\n")

    random_data = {}
    for i in range(5, 51):
        distribution_sum = []
        distribution_avg = []
        for j in range(0,k):
            sample = random.choices(list(interactions.keys()), k=i)
            module_genes = make_artifical_module(sample, i)
            #print(module_genes)
            module_dPCC = []

            for gene1 in range(0, len(module_genes)):
                for gene2 in range(gene1+1, len(module_genes)):
                    g1 = module_genes[gene1]
                    g2 = module_genes[gene2]
                    #print('i')
                   # check if genes have an interaction                       
                    if g1+'_'+g2 in interactions and g1 in genes and g2 in genes:
                            #print('in')
                        g1_exp_case = [genes[g1][t] for t in range(0,len(header)) if header[t] in case]
                        g1_exp_cont = [genes[g1][t] for t in range(0,len(header)) if header[t] in control]
                        g2_exp_case = [genes[g2][t] for t in range(0,len(header)) if header[t] in case]
                        g2_exp_cont = [genes[g2][t] for t in range(0,len(header)) if header[t] in control]


                        PCC_case = np.corrcoef([float(i)+0.000000001 for i in g1_exp_case], [float(i)+0.000000001 for i in g2_exp_case])
                        PCC_cont = np.corrcoef([float(i)+0.000000001 for i in g1_exp_cont], [float(i)+0.000000001 for i in g2_exp_cont])
                        #print(PCC_case[0][1], PCC_cont[0][1])
                        if np. isnan(PCC_case[0][1]) == True and np. isnan(PCC_cont[0][1]) == True:
                            module_dPCC.append(0)
                            dPCC = 0
                        elif np. isnan(PCC_case[0][1]) == True and np. isnan(PCC_cont[0][1]) == False:
                            module_dPCC.append(abs(0 - PCC_cont[0][1]))
                            dPCC = abs(0 - PCC_cont[0][1])
                        elif np. isnan(PCC_case[0][1]) == False and np. isnan(PCC_cont[0][1]) == True:
                            module_dPCC.append(abs(PCC_case[0][1] - 0))   
                            dPCC = abs(PCC_case[0][1] - 0)
                        else:
                            dPCC = abs(PCC_case[0][1] - PCC_cont[0][1])
                            module_dPCC.append(dPCC)
                        
            distribution_sum.append(np.sum(module_dPCC))
            distribution_avg.append(np.mean(module_dPCC))
            random_data[i] = distribution_avg


    # iterate through the modules
    cc= 0
    for module in modules:
        module_genes = modules[module]
        if 'COM' in module:
            continue
        #print(module_genes)
        #exit()
        # continue in case of too large or too small modules
        if len(module_genes) < 5 or len(module_genes) > 50:
            continue


    # iterate through the experiments
        for experiment in experiments:
         
            conditions = experiments[experiment]

            
            # check existence of control and length of control. Skip if not sufficient
            if 'control' in conditions:
                control = conditions['control']
            else:
                continue

            if len(control) < 3:
                continue 

            # iterate throug the conditions   
            for condition in conditions:
                if condition == 'control':
                    continue   
                
                case = conditions[condition]               
                
                if len(case) < 3:
                    continue 

                module_dPCC = []

                for gene1 in range(0, len(module_genes)):
                    for gene2 in range(gene1+1, len(module_genes)):
                        g1 = module_genes[gene1]
                        g2 = module_genes[gene2]
                        # check if genes are in                         
                        if g1+'_'+g2 in interactions and g1 in genes and g2 in genes:
                            #print('in')
                            g1_exp_case = [genes[g1][t] for t in range(0,len(header)) if header[t] in case]
                            g1_exp_cont = [genes[g1][t] for t in range(0,len(header)) if header[t] in control]
                            g2_exp_case = [genes[g2][t] for t in range(0,len(header)) if header[t] in case]
                            g2_exp_cont = [genes[g2][t] for t in range(0,len(header)) if header[t] in control]


                            PCC_case = np.corrcoef([float(i)+0.000000001 for i in g1_exp_case], [float(i)+0.000000001 for i in g2_exp_case])
                            PCC_cont = np.corrcoef([float(i)+0.000000001 for i in g1_exp_cont], [float(i)+0.000000001 for i in g2_exp_cont])

                            if np. isnan(PCC_case[0][1]) == True and np. isnan(PCC_cont[0][1]) == True:
                                module_dPCC.append(0)
                                dPCC = 0
                            elif np. isnan(PCC_case[0][1]) == True and np. isnan(PCC_cont[0][1]) == False:
                                module_dPCC.append(abs(0 - PCC_cont[0][1]))
                                dPCC = abs(0 - PCC_cont[0][1])
                            elif np. isnan(PCC_case[0][1]) == False and np. isnan(PCC_cont[0][1]) == True:
                                module_dPCC.append(abs(PCC_case[0][1] - 0))   
                                dPCC = abs(PCC_case[0][1] - 0)
                            else:
                                dPCC = abs(PCC_case[0][1] - PCC_cont[0][1])
                                module_dPCC.append(dPCC)


                #if np.isnan(module_dPCC) == False:# and np.isnan(random_data[len(module_genes)]) == False:
                z_score = ( np.mean(module_dPCC) - np.mean(random_data[len(module_genes)]) ) / np.std(random_data[len(module_genes)])+0.001
                p_value = 1- stats.norm.cdf(abs(z_score))
                #else:
                 #   z_score = 'nan'
                  #  p_value = 'nan'
               
                #print(module, experiment, condition, len(module_genes), np.sum(module_dPCC), np.mean(module_dPCC), z_score, p_value)
                outfile.write(module+'\t'+str(experiment)+'\t'+condition+'\t'+str(len(module_genes))+'\t'+str(np.mean(module_dPCC))+'\t'+str(z_score)+'\t'+str(p_value)+'\n')
                              

                
def nPCC(outfolder, header, genes, modules, k):
    """
    Calculate nPCC score per module.
    """

    # create Folder to store 
    if os.path.isdir(outfolder+'/nPCC') == False:
        os.mkdir(outfolder+'/nPCC')

    outfile = open(outfolder+'/nPCC/nPCC.csv', 'w')
    outfile.write("module_type\tmodule_no\tmodule_length\tmean_PCC\tfracion_pairs_with_suitable_expression_value\tz_score\tp_values\n")
    # sampling pool
    pool = pooling(modules)
    
    cc= 0
    random_modules = {}
    for i in range(5, 51):
        distribution = []
        for j in range(0,k):
            sample = random.sample(pool,j)
            PCC_values       = []
            counter          = 1

            for gene1 in range(0, len(sample)):
                for gene2 in range(gene1+1, len(sample)):
                    g1 = sample[gene1]
                    g2 = sample[gene2]
                    counter += 1
                    if g1 in genes and g2 in genes:
                        try: 
                            PCC = np.corrcoef([float(i) for i in genes[g1]], [float(i) for i in genes[g2]])
                            #print(PCC)
                            PCC_values.append(abs(PCC[0][1]))
                            #print(gene1, gene2, PCC[0][1])        
                        except:
                            print("Warning")
            if len(PCC_values) > 0:
                distribution.append(np.mean(PCC_values))
            print("Calculate nPCC ", round((cc/(k*45))*100, 2),"%",end="\r")
            cc += 1
        #print(distribution)
        #print(i, np.mean(distribution),  np.std(distribution))
        random_modules[i] = [np.mean(distribution), np.std(distribution)]  


    # iterate through the modules
    for module in modules:
        module_genes     = modules[module]
        PCC_values       = []
        counter          = 1
        if len(module_genes) < 5 or len(module_genes) > 50:
            continue

        for gene1 in range(0, len(module_genes)):
            for gene2 in range(gene1+1, len(module_genes)):
                g1 = module_genes[gene1]
                g2 = module_genes[gene2]
                counter += 1
                if g1 in genes and g2 in genes:
                    try: 
                        PCC = np.corrcoef([float(i) for i in genes[g1]], [float(i) for i in genes[g2]])
                        PCC_values.append(abs(PCC[0][1]))
                    except:
                        print("Warning")
                    
        try:
            mean = np.mean(PCC_values)
        except:
            mean = 'nan'

        # calculate z-score
        z_score  = (mean - random_modules[len(module_genes)][0] ) / (random_modules[len(module_genes)][1]+0.001)
        p_values = 1- stats.norm.cdf(abs(z_score))   
        #print(module, np.mean(PCC_values), round(len(PCC_values)/counter,3), z_score, p_values)
        outfile.write(module.split('_')[0]+'\t'+module.split('_')[1]+'\t'+str(len(module_genes))+'\t'+ str(np.mean(PCC_values))+'\t'+str(round(len(PCC_values)/counter,3))+'\t'+str(z_score)+'\t'+str(p_values)+'\n')
        #exit()                        

        #print(module_genes)      

def calculateActivity(modules, pvalues, outfolder, k):
    """
    Calculate the biological activity of a particular subnetwork per module based on Ideker et al 2002
    """

    # create random sets of genes with certain size
    pool = pooling(modules)
    random_modules = {}
    for i in range(5, 51):
        distribution = []

        for j in range(0,k):
            sample = random.sample(pool,j)
            module_zi        = []
            for gene in sample:
                if gene in pvalues:    
                    pval = float(pvalues[gene])
                    zi   = norm.ppf((1-pval)+0.000000000000001)
                    module_zi.append(zi)
            
                
            if len(module_zi) > 0:        
                za =  (1/(len(module_zi))**0.5)*sum(module_zi)
                #za =  (1/(len(module_zi)))*sum(module_zi)
                distribution.append(za)
        random_modules[i] = [np.mean(distribution), np.std(distribution)]

    #for i in random_modules:
    #    print(i, random_modules[i])


    # create Folder to store 
    if os.path.isdir(outfolder+'/moduleActivity/') == False:
        os.mkdir(outfolder+'/moduleActivity/')

    outfile = open(outfolder+'/moduleActivity/activity.csv', 'w')
    outfile.write("module_type\tmodule_no\tfull_name\tmodule_length\tz_score\tsubnet_score\n")

    cc = 0
    for module in modules:
        module_genes     = modules[module]
        module_zi        = []
        if len(module_genes) <5 or len(module_genes) > 50:
            cc += 1
            continue    

        for gene in module_genes:
            if gene in pvalues:
                pval = float(pvalues[gene])
                zi   = norm.ppf((1-pval)+0.000000000000001)
                module_zi.append(zi)
                
        if len(module_zi) > 0:        
            za =  (1/(len(module_zi))**0.5)*sum(module_zi)
            #za =  (1/(len(module_zi)))*sum(module_zi)
            # correct za
            sa = (za-random_modules[len(module)][0])/random_modules[len(module)][1]
            outfile.write(module.split('_')[0]+'\t'+module.split('_')[1]+'\t'+module+'\t'+str(len(module_genes))+'\t'+str(za)+'\t'+str(sa)+'\n')
        else:
            za =  0    
        print("Calculate activity ", round((cc/len(modules))*100, 2),"%",end="\r")
        cc += 1
        

if __name__ == '__main__':
    """
    Generates a list of preferable non-redundant motifs for a set of directed and a set of undirected edges
    First argument: 
    """    
    # load arugments
    args = arg_parse()

    # create outfolder if not existing
    try:
        if args.out == None:
            print('Provide path to output folder with -out')
            exit()
        os.mkdir(f"{args.out}/")
    except FileExistsError:
        pass

    # nPCC calculation
    if args.nPCC == True:
        print("\nPreparing nPCC calculation\n")
        
        # load modules
        print("Load modules ",end="\r")
        try:
            modules = loadModules(args.clusters)
            print("Load modules \t\t... Done!")
        except:
            print('Error: Provide correct path to SUBATOMIC SCHYPE folder with -clusters')
        
        # load expression data 
        print("Load expression",end="\r")
        try:
            header, genes = loadExpression(args.expression)
            print("Load expression \t... Done!")
        except:
            print('Error: Provide correct path to exrpession data with -expression')

        # Check parameter k  
        print("Check parameters",end="\r")
        try:
            k = int(args.sampling_k)
            print("Check parameters \t... Done!")
        except:
            print('Error: Parameter -sampling_k not correctly set')

        # Calculation of nPCC
        print("Calculate nPCC",end="\r")
        try:
            nPCC(args.out, header, genes, modules, k)
            print("                                                    ", end = '\r')
            print("Calculate nPCC \t\t... Done!")
        except:
            print('Error: calculation of nPCC failed. Pelase check all input files')       

    # activity calculation
    if args.activity == True:
        print("\nPreparing module activity calculation\n")

        # load modules
        print("Load modules ",end="\r")
        try:
            modules = loadModules(args.clusters)
            print("Load modules \t\t... Done!")
        except:
            print('Error: Provide correct path to SUBATOMIC SCHYPE folder with -clusters')

        # load p-values
        print("Load p-values ",end="\r")
        try:
           pvalues = loadPvalues(args.pvalues)
           print("Load p-values \t\t... Done!")
        except:
            print('Error: Provide correct path to p-value file with -pvalues')

        # Check parameter k  
        print("Check parameters",end="\r")
        try:
            k = int(args.sampling_k)
            print("Check parameters \t... Done!")
        except:
            print('Error: Parameter -sampling_k not correctly set')

        # calculate biological activity
        print("Calculate activity ",end="\r")
        try:
           calculateActivity(modules, pvalues, args.out, k)
           print("                                                    ", end = '\r')
           print("Calculate activity \t... Done!")
        except:
            print('Error: calculation of module activity failed. Pelase check all input files')

    # calculate ECD score
    if args.ECD == True:
        print("\nPreparing ECD calculation\n")    
        
        # load modules
        print("Load modules ",end="\r")
        try:
            modules = loadModules(args.clusters)
            print("Load modules \t\t... Done!")
        except:
            print('Error: Provide correct path to SUBATOMIC SCHYPE folder with -clusters')

        # load modules
        print("Load interacions ",end="\r")
        try:
            interactions = loadInteractions(args.clusters)
            print("Load interactions \t... Done!")
        except:
            print('Error: Provide correct path to SUBATOMIC SCHYPE folder with -clusters')

        # load expression data 
        print("Load expression",end="\r")
        try:
            header, genes = loadExpression(args.expression)
            print("Load expression \t... Done!")
        except:
            print('Error: Provide correct path to exrpession data with -expression')

        # Check parameter k  
        print("Check parameters",end="\r")
        try:
            k = int(args.sampling_k)
            print("Check parameters \t... Done!")
        except:
            print('Error: Parameter -sampling_k not correctly set')

        # Load meta file  
        print("Load meta data",end="\r")
        try:
            experiments, case, control = loadMeta(args.meta)
            print("Load meta data \t\t... Done!")
        except:
            print('Error: path to meta data or meta data syntax not correctly set')

        # calculate biological activity
        print("Calculate ECD ",end="\r")
        try:
            ECDscore(args.out, experiments, case, control, header, genes, modules, interactions, k)
            print("                                                    ", end = '\r')
            print("Calculate ECD \t\t... Done!")
        except:
            print('Error: calculation of module activity failed. Pelase check all input files')

    print('\nAll calcualtions are completed\n')



     
    

    
