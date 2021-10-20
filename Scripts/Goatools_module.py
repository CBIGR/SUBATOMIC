"""
@author: hayoung kim
"""

from sys import argv
import subprocess
import os


def make_dir(path):
    """
    Makes new directories
    """
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


def parseArguments(argvs):
    """
    This function takes the comandline arguments of argv and returns
    a dictionary that contains these arguments with flag as key
    """
    parameter = {}

    for i in range(1, len(argvs)):
        # load motiftype
        if argvs[i] == "--type" or argvs[i] == "-t":
            parameter["type"] = argvs[i + 1]
        # load file that contains the motifs to search for
        if argvs[i] == "--modules" or argvs[i] == "-m":
            parameter["modules"] = argvs[i + 1]
        # folder to write output
        if argvs[i] == "--out" or argvs[i] == "-o":
            parameter["out"] = argvs[i + 1]
        if argvs[i] == "--species" or argvs[i] == "-s":
            parameter["species"] = argvs[i + 1]
        # relative or absolute path to the file containing the GO associations
        if argvs[i] == "--go" or argvs[i] == "-g":
            parameter["GO_associations"] = argvs[i + 1]
        # current run file of the parallel split
        if argvs[i] == "--run" or argvs[i] == "-r":
            parameter["run"] = argvs[i + 1]
        # path to annotation file or mode of annotation (ModuleType, All, Custom)
        if argvs[i] == "--annotation" or argvs[i] == "-a":
            parameter["annotation"] = argvs[i + 1]
        # path to the full network interaction file
        if argvs[i] == "--interactions" or argvs[i] == "-i":
            parameter["interactions"] = argvs[i + 1]
        if argvs[i] == "--pvalue" or argvs[i] == "-p":
            parameter["pvalue"] = argvs[i + 1]
    return parameter


def loadModuleDic(path):
    """
    load modules and store them in a dictionary
    """
    modules_dict = {}
    infile = open(path, 'r')
    lines = infile.read().split('\n')

    for line in lines:
        tmp = line.split('\t')
        if len(tmp) > 0 and len(line) > 0:
            key = tmp[1]
            modules_dict.setdefault(key, [])
            modules_dict[key].append(tmp[0])
    return modules_dict


def MakeModuleTextFile(species, motif_type, module):
    """
    Make text files classified by each module
    """

    with open(f'{species}/goatools-main/data/module_textfile/SCHype{motif_type}/module_{str(module)}.txt', 'w+') as module_text_file:
        for gene in Modules[module]:
            genes = gene + '\n'
            module_text_file.write(genes)


def new_nodes_file(params, motif_type, rID):
    """
    Makes the right input format from the nodes file
    """
    
    species = params["species"]
    output = open(f"{species}/goatools-main/data/new_{motif_type}_{rID}.nodes.txt", "w+")
    with open(f"{species}/SCHYPE/SCHype{motif_type}/{motif_type}.nodes.txt") as i:
        for line in i:            
            line = line.rstrip()
            accession, cluster = line.split('\t')
            output.write(f"{accession}\n")

    output.close()

def new_nodes_file2(params, motif_type, rID):
    """
    Makes the right input format from whole network file
    """
    node_dict = {}
    species = params["species"]
    output = open(f"{species}/goatools-main/data/new_{motif_type}_{rID}.nodes.txt", "w+")
    with open(param["interactions"]) as i:
        for line in i:            
            line = line.rstrip()
            tmp = line.split('\t')
            if len(tmp) > 2:
                node_dict[tmp[0]] = ''
                node_dict[tmp[1]] = ''
    
    for i in node_dict:
        output.write(f"{i}\n")

    output.close()

def new_nodes_file3(params, motif_type, rID):
    """
    Load directly from file
    """
    species = params["species"]
    output = open(f"{species}/goatools-main/data/new_{motif_type}_{rID}.nodes.txt", "w+")
    with open(params["annotation"]) as i:
        for line in i:            
            line = line.rstrip()
            output.write(f"{line}\n")

    output.close()



def runGOATOOLS(params, motif_type, modules, go_associations, rID, loadType, p_value):
    """
    Run GOATOOLS per every module
    """
    # counter for motifs
    r = 1
    species = params["species"]
    # go through every single motifs:
    for module in modules.keys():
        MakeModuleTextFile(species, motif_type, module)
        
        make_dir(f"{species}/GOATOOLS/SCHype{motif_type}")
       

        goatools_run = f'python3 Scripts/Goatools_find_enrichment.py {species}/goatools-main/data/module_textfile/SCHype{motif_type}/module_{str(module)}.txt {species}/goatools-main/data/new_{loadType}_{rID}.nodes.txt {go_associations} --pval='+str(p_value)+f' --method=fdr_bh --pval_field=fdr_bh --outfile={params["species"]}/GOATOOLS/SCHype{motif_type}/{str(module)}.txt'


        # execution of command line
        subprocess.run(goatools_run, shell=True)
        r += 1

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

def filterModules(modules, filterSet):
    """
    Only keep modules that are designated for this parallel run
    """
    new_modules = {}

    for module in modules:
        if module in filterSet:
            new_modules[module] = modules[module]
    return new_modules
    

if __name__ == '__main__':
    """
    This script defines motiftype and
    processes outputs module dictionary
    """
    # load all command line arguments and store them in the param dict
    param = parseArguments(argv)
    # load motif type
    Motif_Type = param["type"]
    # load split
    split = readParallelSplit(param["run"])
    # get run identifier
    rID = param["run"].replace('\\',"/").split('/')[-1].split(".")[0]


    # make the necessaery folders
    for i in split.keys():

        Species = param["species"]
        make_dir(f"{Species}")
        make_dir(f"{Species}/goatools-main")
        make_dir(f"{Species}/GOATOOLS")
        make_dir(f"{Species}/goatools-main/data")
        make_dir(f"{Species}/goatools-main/data/module_textfile")
        make_dir(f"{Species}/goatools-main/data/module_textfile/SCHype{i}")
    
        # make the input for the nodes file
        if param["annotation"] == "ModuleType": 
            new_nodes_file(param, i, rID)
        elif param["annotation"] == "All":
            new_nodes_file2(param, "full", rID)
        else:
            new_nodes_file3(param, "full", rID)


        # load module dictionary
        Modules = loadModuleDic(f'{param["modules"]}/SCHype{i}/{i}.nodes.txt')
                
        # Filter modules
        Modules = filterModules(Modules, split[i])

        # get background
        if param["annotation"] == "ModuleType": 
            runGOATOOLS(param, i, Modules, param["GO_associations"], rID, i, float(param["pvalue"]))
        else:
            runGOATOOLS(param, i, Modules, param["GO_associations"], rID, "full", float(param["pvalue"]))


    try:
        os.mkdir(param["species"]+"/Logs")
    except FileExistsError:
        pass
    with open(param["species"]+"/Logs/Enrichment_"+param["species"]+"_"+rID+"_done.txt", 'w+') as f:
        f.write("Enrichment view done!")
# CLI
# python goatools_module.py -t COM -m data/ -o Atha_goatools_outputs/SCHypeCOM
