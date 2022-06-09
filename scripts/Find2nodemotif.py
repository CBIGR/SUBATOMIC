# Find 2-node motifs
#######################
import argparse
import os


def arg_parse():
    """
    Parses the arguments given in the command line
    """
    parser = argparse.ArgumentParser(description='Find 2node motifs based on interactionfiles.')
    parser.add_argument('input', help='The folder where your input files are located')
    parser.add_argument('output', help='The folder where your output files need to be written to')
    parser.add_argument('directed_interactions', help='List of directed interactions')
    parser.add_argument('undirected_interactions', help='List of undirected interactions')

    return parser.parse_args()


def undirected(folder, interaction):
    """
    Opens all undirected interactions and puts them in a specified list.
    """
    with open(f'{folder}/{interaction}', 'r') as f:
        for line in f:
            line = line.rstrip()
            tmp = line.split('\t')
            U.append(f'{tmp[0]}\t{tmp[1]}\t{interaction}')
            U.append(f'{tmp[1]}\t{tmp[0]}\t{interaction}')


def directed(folder, interaction):
    """
    Opens all directed interactions and puts them in a specified list.
    """
    with open(f'{folder}/{interaction}', 'r') as f:
        for line in f:
            line = line.rstrip()
            tmp = line.split('\t')
            D.append(f'{tmp[0]}\t{tmp[1]}\t{interaction}')


def find_2_node_du_motif(list_directed, list_undirected, out_folder):
    """
    Finds all 2-node DU motifs.
    The file "DU2nodemotifs.txt" is automatically created in the given output folder.
    If it already exists, the new info will be appended to the existing file.
    """
    with open(f'{out_folder}/DU2nodemotifs.txt', 'w') as f:
        for d in list_directed:
            col_d = d.split('\t')
            for u in list_undirected:
                col_u = u.split('\t')
                if col_d[0] == col_u[0] and col_d[1] == col_u[1]:
                    f.write(f'{col_d[0]}\t{col_d[1]}\t{col_d[2][-7]}{col_u[2][-7]}\n')


def find_2_node_dd_motif(list_directed, out_folder):
    """
    Finds all 2-node DD motifs.
    The file "DD2nodemotifs.txt" is automatically created in the given output folder.
    If it already exists, the new info will be appended to the existing file.
    """
    with open(f'{out_folder}/DD2nodemotifs.txt', 'w') as f:
        for d in list_directed:
            col_d = d.split('\t')
            for d2 in list_directed:
                col_d2 = d2.split('\t')
                if col_d[1] == col_d2[0] and col_d[0] == col_d2[1]:
                    f.write(f'{col_d[0]}\t{col_d[1]}\t{col_d[2][-7]}{col_d2[2][-7]}\n')


def open_input(infolder, undirected_interactions, directed_interactions):
    """
    Opens all the files in the given input folder and filters them based on the fact if the interaction files are
    directed or undirected.
    The files will be opened by their specific function.
    """
    ui = undirected_interactions.split('_')                   # Contains all undirected interactions
    di = directed_interactions.split('_')      # Contains all directed interactions
    for f in os.listdir(infolder):
        if f[-7] in ui:
            undirected(infolder, f)
        elif f[-7] in di:
            directed(infolder, f)


def make_dir(name):
    """
    Makes new directories.
    If the directory already exists, it escapes that error.
    """

    try:
        os.mkdir(name)
    except FileExistsError:
        pass


if __name__ == '__main__':
    """
    Finds the 2 node motifs, either DD or DU (specified in the command line).
    The input files are the interaction files.
    The output file is either the DD2nodemotif.txt or the DU2nodemotif.txt .
    """

    # Parse the command line arguments
    args = arg_parse()

    # Check if the given input path exists
    assert os.path.isdir(args.input), f'{args.input} not found'

    # The input is opened and the directed and undirected interactions are put in their respective lists
    U = []
    D = []
    open_input(args.input, args.undirected_interactions, args.directed_interactions)

    # Check if the given output path exists
    for i, file in enumerate(args.output.split('/')):
        if i == 0:
            make_dir(file)
        else:
            make_dir(f'{"/".join(args.output.split("/")[:i])}/{file}')
    assert os.path.isdir(args.output), f'{args.output} not found'

    # The DU or the DD subgraphs are made
    find_2_node_du_motif(D, U, args.output)
    print("DU subgraphs made")

    find_2_node_dd_motif(D, args.output)
    print("DD subgraphs made")
