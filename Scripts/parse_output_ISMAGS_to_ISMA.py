"""
@author: jens loers
"""

from sys import argv
import subprocess
from os import listdir, mkdir, system
from os.path import isfile, join


def parseArguments(argv):
    """
    This function takes the comandline arguments of argv and returns
    a dictionary that contains these arguments with flag as key
    """
    parameter = {}

    for i in range(1, len(argv)):
        # path to ISMAGS output
        if argv[i] == "--folder" or argv[i] == "-f":
            parameter["folder"] = argv[i + 1]
        if argv[i] == "--species" or argv[i] == "-s":
            parameter["species"] = argv[i + 1]

    return parameter

def getFilenames(path):
    """
    Returns all file names from a given folder
    """
    return [f for f in listdir(path) if isfile(join(path, f))]

def parser(path):
    """
    Conversion of ISMAGS output into ISMA output
    """
    
    files = getFilenames(path)

    for f in files:
        # get motif name
        motif = f[:3]
        # store lines
        infile = open(path+'/'+f, 'r')
        lines  = infile.read().split('\n')
        infile.close()
        new_lines = []
        for line in lines:
            tmp = line.split(';')
            if len(tmp) > 1:
                new_lines.append("Motif ["+motif+"]: "+"["+tmp[0][:-1]+", "+tmp[1][:-1]+", "+tmp[2][:-1]+"]")

        outfile = open(path+'/'+f, 'w')
        outfile.write(str(len(new_lines))+' motifs found\n')
        for line in new_lines:
            outfile.write(line+'\n')

#        print(motif)


if __name__ == '__main__':
    """
    This script formats the ISMAGS out into the ISMA output and removes
    the letters ISMA adds to the 
    """

    # load all command line arguments and store them in the param dict
    param   = parseArguments(argv)
    species = param["folder"] 
    # start parsing
    parser(param["species"])

    with open(f'{param["species"]}/Logs/ISMAGS_ISMA_PARSE_done.txt', 'w+') as f:
        f.write("ISMAGS_ISMA_PARSE_done done!")
