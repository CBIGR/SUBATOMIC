"""
Created on 2020. 11. 2.
@author: jens loers
"""

from sys import argv
import subprocess
from os import listdir, mkdir, system
from os.path import isfile, join


def make_dir(name):
    """
    Makes new directories.
    If the directory already exists, it escapes that error.
    """

    try:
        mkdir(name)
    except FileExistsError:
        pass


def parseArguments(argv):
    """
    This function takes the comandline arguments of argv and returns
    a dictionary that contains these arguments with flag as key
    """
    parameter = {}

    for i in range(1, len(argv)):
        # path to ISMA
        if argv[i] == "--ISMAGS" or argv[i] == "-i":
            parameter["ISMAGS"] = argv[i + 1]
        # load file that contains the motifs to search for
        if argv[i] == "--motifs" or argv[i] == "-m":
            parameter["motifs"] = argv[i + 1]
        # path to folder containing ONLY the interaction files
        if argv[i] == "--interaction_folder" or argv[i] == "-f":
            parameter["folder"] = argv[i + 1]
        # folder to write output
        if argv[i] == "--out" or argv[i] == "-o":
            parameter["out"] = argv[i + 1]
        # the organism you are working in
        if argv[i] == "--species" or argv[i] == "-s":
            parameter["species"] = argv[i + 1]

    return parameter


def loadMotifFIle(path):
    """
    load motifs and store them in a list
    """
    motifs = []
    infile = open(path, 'r')
    lines = infile.read().split('\n')

    for line in lines:
        tmp = line.split('\t')
        if len(tmp) > 0 and len(line) > 0:
            motifs.append(tmp[0])
    return motifs


def makeInteractionFileSyntax(folder):
    """
    Careful: Files need to be in one folder with ONLY the interaction files as input
    NAMING convetion: Filename_interactionName_directed.txt
    EXAMPLE: Atha_P_u.txt
    Only with that information in the file name, the input can be dynamic
    """

    files = {}

    # get all files names in the input folder
    onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]

    # create the right syntax for ISMAs linkfile input
    for name in onlyfiles:
        #sp, interaction, direction 
        names        = name.split("_")
        sp           = names[-1]
        interaction  = names[-2]
        direction    = names[-3]
        direction    = direction.split(".")[0]
        files[interaction.upper()] = interaction.upper() + " " + direction + " a a " + name
        files[interaction.lower()] = interaction.lower() + " " + direction + " a a " + name

    return files


def runISMA(param, motifs, syntax):
    """
    Run ISMA per motif
    """
    # counter for motifs
    R = 1

    # go through every single motifs:
    for motif in motifs:
        # link_files are the files that need to be loaded to calculate the motifs
        #print(syntax)
        #link_files = '"' + syntax[motif[0]] + " " + syntax[motif[1]] + " " + syntax[motif[2]] + '"'

        # select only files for motifs
        used = []
        link_files = '"'
        for m in motif:
            if m not in used:
                link_files += syntax[m] +  " "
                used.append(m)
        link_files = link_files[:-1] + '"'
                     

        # command line call to run ISMA
        isma_run = f'java ' \
                   f'-jar {param["ISMAGS"]} ' \
                   f'-folder {param["folder"]}/ ' \
                   f'-linkfiles {link_files} ' \
                   f'-motif {motif} ' \
                   f'-output {param["out"]}/{motif}-{R}.txt'

        # execution of command line
        #print(isma_run)
        #system(isma_run)
        R += 1
    #exit()

if __name__ == '__main__':
    """
    This script processes interaction files and outputs motif files
    for a list of pre-defied motifs 
    """

    # load all command line arguments and store them in the param dict
    param = parseArguments(argv)

    make_dir(f'{param["species"]}/Logs')
    make_dir(f'{param["species"]}/3nodemotifs')


    # load motif files
    motifs = loadMotifFIle(param["motifs"])
    # make interaction file sythax
    #print(param["folder"])
    syntax = makeInteractionFileSyntax(param["folder"])
    # run ISMA
    runISMA(param, motifs, syntax)

    with open(f'{param["species"]}/Logs/ISMA_done.txt', 'w+') as f:
        f.write("ISMAGS done!")
