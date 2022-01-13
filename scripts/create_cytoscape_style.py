import sys
import os
import csv
import pandas as pd
import itertools
import os
import sys
import argparse


def arg_parse():
    """
    Parses the arguments given in the command line
    """
    parser = argparse.ArgumentParser(description='Collect results from parallelization.')
    parser.add_argument('-folder', help='The standard run folder.')
    parser.add_argument('-cytoscape', help='The relative or absolute path to folder where the standard_cytoscape style sheet is located.')
    parser.add_argument('-di', help='Directed interactions that need an arrow')
    parser.add_argument('-ui', help='Directed interactions that need an arrow')
    return parser.parse_args()
    
def loadStyle(path):
    """
    Load the Cytoscape style xml
    """

    infile = open(path, 'r')
    return infile.read().split('\n')    
    

if __name__ == '__main__':
    """
    Modify the cytoscape style sheet to know which network letters are directed
    """

    #parse argumnets
    args = arg_parse()    
    # get directed interaction list
    di = args.di.split('_')
    # get undirected interaction list
    ui = args.ui.split('_')
    # load stype
    style =  loadStyle(args.cytoscape)
    # write new stype file
    if not os.path.exists('{args.folder}/Cytoscape/'):
        os.makedirs(args.folder+'/Cytoscape')

    outfile = open(args.folder + '/Cytoscape/cytoscape_style.xml', 'w')
    skip = False
    for i in range(0, len(style)):
        if skip == True:
            skip = False
            continue
        if style[i] == '            <visualProperty default="NONE" name="EDGE_TARGET_ARROW_SHAPE">':
            outfile.write(style[i]+'\n')
            outfile.write(style[i+1]+'\n')
            i+= 2
            skip = True
            for j in di:
                 outfile.write('                    <discreteMappingEntry attributeValue="'+j+'" value="DELTA"/>'+'\n')
            for j in ui:
                 outfile.write('                    <discreteMappingEntry attributeValue="'+j+'" value="NONE"/>'+'\n')
        else:
            outfile.write(style[i]+'\n')



    try:
        os.mkdir(f"{args.folder}/Logs")
    except FileExistsError:
        pass
    with open(f"{args.folder}/Logs/Cytoscape_style_done.txt", 'w+') as f:
        f.write("Enrichment collector done!")
