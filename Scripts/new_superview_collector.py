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
    parser = argparse.ArgumentParser(description='Collect results from parallelization.')
    parser.add_argument('-folder', help='The relative or absolute path to folder where your input files are located, in this folder, also the output will be written.')
    parser.add_argument('-runs', help='Current run of parallelization')
    parser.add_argument('-groups', help='Modules bzw groups')
    #parser.add_argument('Path_to_SCHYPE_jar', help='The relative or absolute path to the master-java folder containing the SCHYPE java script.')

    return parser.parse_args()


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
    runs = args.runs.split('_')

    #print(runs)
    

    # groups to run against
    groups = args.groups.split('_')
    
    # collect per run
    header = ''
        
    store = {}

    for run in runs:
        path = args.folder+'/Superview/'+run+'.txt'

         
        ## calculate superview
        #name = run.replace("\\","/").split('/')[-1]
        #path = args.folder+'/Superview/'+name
        try:
            infile = open(path, 'r')
        except:
            continue
        lines  = infile.read().split('\n')
        #print(lines)
        header = lines[0]
    
        for line in lines[1:]:
            tmp = line.split('\t')
            if len(line) > 1:
                if tmp[0] not in store:
                    store[tmp[0]] = [line]
                else:
                    store[tmp[0]].append(line)
        #print(store)
                               
    for key in store:
        outfile = open(args.folder+'/Superview/'+key+'.csv', 'w')
        outfile.write(header+'\n')
        for line in store[key]:
             outfile.write(line+'\n')
    
     

    for run in runs:
        path = args.folder+'/Superview/'+run+'.txt'
        if os.path.exists(path):
            os.remove(path)




try:
    os.mkdir(args.folder+"/Logs")
except FileExistsError:
    pass
with open(args.folder+"/Logs/Super_View_Collector_done.txt", 'w+') as f:
    f.write("Super view done!")



    
