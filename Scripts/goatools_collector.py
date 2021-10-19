import sys
import math
import os
import sys
import argparse
from os import listdir
from os.path import isfile, join

def arg_parse():
    """
    Parses the arguments given in the command line
    """
    parser = argparse.ArgumentParser(description='Collect results from parallelization.')
    parser.add_argument('-folder', help='The relative or absolute path to folder where your input files are located, in this folder, also the output will be written.')
    parser.add_argument('-runs', help='Current run of parallelization')
    parser.add_argument('-groups', help='Modules bzw groups')
    #parser.add_argument('Path_to_SCHYPE_jar', help='The relative or absolute path to the master-java folder containing the SCHYPE java script.')

    return parser.parse_args()

def rankPvalues(pvalues):
    """
    return for a dinctionary with p-values their increasing rank
    """
    rank = 1
    
    for i in sorted(pvalues.keys()):
        pvalues[i] = rank
        rank += 1

    return pvalues


if __name__ == '__main__':
    """
    Script to collect all outputs of the Go-annotations into one file
    """

    #parse argumnets
    args = arg_parse()
    

    # load split for parallelization
    runs = args.runs.split('_')


    # groups to run against
    groups = args.groups.split('_')
    
    # collect per run
    #header = ''
        
    store = {}
    header = ''

    content_storage = []

    for group in groups:
        #/home/jloers/repo/composite_subgraph_cluster_integration/Celegans/GOATOOLS/SCHypeCIR
        path = args.folder+'/GOATOOLS/SCHype'+group
        #print(path, os.path.isfile(path))
        #exit()
        #if os.path.isfile(path) == False:
        #    continue
         
        ## calculate superview
        #name = run.replace("\\","/").split('/')[-1]
        #path = args.folder+'/Superview/'+name
        
        #infile = open(path, 'r')
        #lines  = infile.read().split('\n')
        #print(lines)
        

       
        try:
            onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        except:
            continue
        
        #print(onlyfiles)
        #exit()
        # for file in folder:
        #print(onlyfiles)

        # read content from every single file and store it properly in a list
        
        for f in onlyfiles:
            infile = open(path+'/'+f, 'r')
            rID    = f.split('.')[0]
            lines  = infile.read().split('\n')
            header = ["group", "module_no", "p-value-rank", "log2foldRatio"]+lines[0].split('\t')
            
            # get p-value ranks
            pvalues = {}
            for line in lines[1:]:
                tmp = line.split('\t')
                if len(tmp) > 1:
                    pvalues[float(tmp[9])] = 0
            pvalues = rankPvalues(pvalues)            


            for line in lines[1:]:
                tmp      = line.split('\t')
                if len(tmp) > 1:
                    try:
                        log2fold = round(math.log((float(tmp[4].split('/')[0])/float(tmp[4].split('/')[1]))/(float(tmp[5].split('/')[0])/float(tmp[5].split('/')[1]))+1,2),3)
                        print(log2fold)
                    except:
                        log2fold = NA
                    content_storage.append([group, rID, str(pvalues[float(tmp[9])]), str(log2fold)]+tmp)

    # Create a file containing all enrichments and calculate enrichment        
    outfile =  open(args.folder+'/GOATOOLS/Enrichment_Modules.csv', 'w')
    outfile.write("\t".join(header)+'\n')
    for entry in content_storage:
        outfile.write("\t".join(entry)+'\n')

    
    

try:
    os.mkdir(args.folder+"/Logs")
except FileExistsError:
    pass
with open(args.folder+"/Logs/Enrichment_Collector_done.txt", 'w+') as f:
    f.write("Enrichment collector done!")



    
