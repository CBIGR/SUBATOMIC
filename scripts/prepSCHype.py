import os
import re
import shutil
import argparse


def arg_parse():
    """
    Parses the arguments given in the command line
    """
    parser = argparse.ArgumentParser(description='Prep and execute SCHYPE.')
    parser.add_argument('folder', help='The relative or absolute path to folder where your input files are located, in this folder, also the output will be written.')
    parser.add_argument('P_value_SCHYPE', default=1, type=int, help='Value calculates ratio between edges and nodes.')
    parser.add_argument('directed_interactions', help='List the directed interaction letters.')
    parser.add_argument('undirected_interactions', help='List the undirected interaction letters.')
    parser.add_argument('Path_to_SCHYPE_jar', help='The relative or absolute path to the master-java folder containing the SCHYPE java script.')

    return parser.parse_args()


def remove_suffix(filename):
    """
    removes suffixes on the accession numbers created by ISMA
    """
    with open(filename, 'r') as i:
        lines = i.readlines()
    with open(filename, 'w') as o:
        for rule in lines:
            new_line = rule.replace('_1', '')
            o.write(new_line)


def found_3_node_motif_file(path, files, directed_interactions, undirected_interactions):
    """
    Writes the found3nodemotifNR.txt file, based on the subgraphs made by ISMA
    Found3nodemotifsNR.txt contains on each line:
        motif-nr.txt \t subgraph \t type of the subgraph
    """

    with open(f"{path}/Found3nodemotifsNR.txt", "w") as o:
        for file_name in files:
            motif_name = file_name[:3]

            # map all diretced and undirected edges to placeholder letters
            network_dict = {}
            for i in directed_interactions:
                network_dict[i.upper()] = 'D'
                network_dict[i.lower()] = 'd'
            for i in undirected_interactions:
                network_dict[i] = 'B'

            new_name = motif_name.translate(
                str.maketrans(network_dict))
            type1 = re.sub('DDD|DDd|dDD|Ddd|ddD|ddd', 'FFL', new_name)
            type2 = re.sub('DdD|dDd', 'CIR', type1)
            type3 = re.sub('BDD|DBd|ddB', 'COP', type2)
            type4 = re.sub('BDd|BdD|DBD|dBd|DdB|dDB', 'FBU', type3)
            type5 = re.sub('Bdd|dBD|DDB', 'COR', type4)
            type6 = re.sub('BBD|BBd|BDB|BdB|DBB|dBB', 'FB2U', type5)
            type7 = re.sub('BBB', 'COM', type6)
            o.write(f"{file_name}\t{motif_name}\t{type7}\n")


def make_type_hash(folder, schype):
    """
    Collects the type of subgraphs found via ISMA
    Contains the files from ISMA per group
    """
    type_hash = {}
    with open(f"{folder}/Found3nodemotifsNR.txt", "r") as motifs:
        for file in motifs:
            file = file.rstrip()
            if file.split('\t')[2] not in type_hash:
                type_hash[file.split('\t')[2]] = []
            type_hash[file.split('\t')[2]].append(file)

    return type_hash


def motifs_type_all(folder, schype):
    """
    Collects the subgraph numbers and motifs generated by ISMA for all subgraphs (ALL group)
    Writes this to the MotifsType.txt file for the ALL group
    """
    motif_re = re.compile(r"\[(.*?)\]")         # regular expression to find the subgraphs in the ISMA output
    with open(f"{folder}/Found3nodemotifsNR.txt", "r") as motifs:
        for file in motifs:
            file = file.rstrip()
            motif_nr, motif = file.split("\t")[:2]
            with open(f"{folder}/3nodemotifs/{motif_nr}", "r") as input_all:            # open subgraph file from ISMA
                with open(f"{folder}/{schype}/SCHypeALL/MotifsType.txt", "a") as output_all:
                    for rule in input_all:
                        rule = rule.rstrip()
                        for match in motif_re.finditer(rule):
                            if len(match.group()) > 5:
                                sub_graph = match.group()[1:-1].split(', ')
                                output_all.write(f"{sub_graph[0]}\t{sub_graph[1]}\t{sub_graph[2]}\t{motif}\n")


def motif_file(folder, type_hash, schype):
    """
    Makes subgraph files for all the different types, with the exception of the ALL group
    Files created:
        MotifsType.txt:     node1 \t node2 \t node3 \t subgraph
        Motifs.txt:         node1 \t node2 \t node3
    """
    motif_re = re.compile(r"\[(.*?)\]")          # regular expression to find the subgraphs in the ISMA output
    for type_folder in type_hash:
        motifs = type_hash[type_folder]          # type hash contains the respective files for the group

        # make a folder for the group
        try:
            os.mkdir(f"{folder}/{schype}/SCHype{type_folder}")
        except FileExistsError:
            pass

        for file in motifs:
            file = file.rstrip()
            motif_nr, motif = file.split("\t")[:2]

            with open(f"{folder}/3nodemotifs/{motif_nr}", "r") as input_group:
                # both files are created at the same moment
                output = open(f"{folder}/{schype}/SCHype{type_folder}/MotifsType.txt", "a")
                output2 = open(f"{folder}/{schype}/SCHype{type_folder}/Motifs.txt", "a")
                for rule in input_group:
                    rule = rule.rstrip()
                    for match in motif_re.finditer(rule):
                        if len(match.group()) > 5:
                            sub_graph = match.group()[1:-1].split(', ')
                            output.write(f'{sub_graph[0]}\t{sub_graph[1]}\t{sub_graph[2]}\t{motif}\n')
                            output2.write(f'{sub_graph[0]}\t{sub_graph[1]}\t{sub_graph[2]}\n')
                output.close()
                output2.close()


def execute(folder, schype, motif_type, path_to_jar, p_value_schype):
    """
    Executes the SCHype java script for the different groups
    """
    # The ALL group has another input file, which also contains the 2-node subgraphs
    if motif_type == 'ALL':
        hg_file = f"{args.folder}/{SCHYPE}/SCHypeALL/MotifsALL.txt"
    else:
        hg_file = f"{folder}/{schype}/SCHype{motif_type}/Motifs.txt"

    shell_command = f'java ' \
                    f'-jar {path_to_jar} ' \
                    f'-hgfile {hg_file} ' \
                    f'-output {folder}/{schype}/SCHype{motif_type}/{motif_type} ' \
                    f'-p {p_value_schype}\n'
    print(f'\n\n{shell_command}\n\n')
    os.system(shell_command)


def execute_momo(folder, schype, path_to_jar):
    """
    Executes the SCHype java script on every single subgraph generated by ISMA
    """
    # Make the MOMO folder
    try:
        os.mkdir(f"{folder}/{schype}/SCHypeMOMO")
    except FileExistsError:
        pass

    # Make a separate folder to make the correct layout on the input folders for SCHype
    try:
        os.mkdir(f"{folder}/{schype}/SCHypeMOMO/subgraphs")
    except FileExistsError:
        pass

    motif_re = re.compile(r"\[(.*?)\]")     # regular expression for finding the motifs in the ISMA output
    motif_files = [x for x in os.listdir(f"{folder}/3nodemotifs")]  # list of all the files coming from ISMA

    with open(f"{folder}/{schype}/SCHypeMOMO/Motiffiles.txt", "w") as motiffile:
        for motif in motif_files:
            # make the motiffiles.txt file (on each line the name of the file in the ISMA folder)
            motiffile.write(motif + '\n')

            # make the correct input format file for SCHype
            with open(f"{folder}/3nodemotifs/{motif}", "r") as i:
                with open(f"{folder}/{schype}/SCHypeMOMO/subgraphs/{motif}", "w") as o:
                    for rule in i:
                        rule = rule.rstrip()
                        for match in motif_re.finditer(rule):
                            if len(match.group()) > 5:
                                sub_graph = match.group()[1:-1].split(', ')
                                o.write(f"{sub_graph[0]}\t{sub_graph[1]}\t{sub_graph[2]}\n")

            hg_file = f"{folder}/{schype}/SCHypeMOMO/subgraphs/{motif}"     # input file for SCHype

            motif = re.sub(r"\.txt", "", motif)
            motif = re.sub("-", "", motif)

            shell_command = f"java -jar {path_to_jar} -hgfile {hg_file} -output {folder}/{schype}/SCHypeMOMO/{motif}\n"
            os.system(shell_command)

    shutil.rmtree(f"{folder}/{schype}/SCHypeMOMO/subgraphs")


if __name__ == '__main__':

    # parse the command line arguments
    args = arg_parse()

    # create the correct version of the SCHype
    SCHYPE = 'SCHYPE'
    #if args.P_value_SCHYPE != 1:
    #    SCHYPE += str(args.P_value_SCHYPE)

    # make folder for the SCHype output
    try:
        shutil.rmtree(f"{args.folder}/{SCHYPE}")
    except:
        pass
    try:
        os.mkdir(f"{args.folder}/{SCHYPE}")
    except FileExistsError:
        pass

    # remove the suffixes on the accession numbers created by ISMA
    #for File in os.listdir(f"{args.folder}/3nodemotifs"):
    #    remove_suffix(f"{args.folder}/3nodemotifs/{File}")

    # create the Foun3nodemotifsNR.txt file
    found_3_node_motif_file(args.folder, os.listdir(f"{args.folder}/3nodemotifs"), args.directed_interactions.split('_'), args.undirected_interactions.split('_') )
    
    # make the folder for the ALL group
    try:
        os.mkdir(f"{args.folder}/{SCHYPE}/SCHypeALL")
    except FileExistsError:
        pass

    # make the type hash containing the files corresponding to each group
    Type_hash = make_type_hash(args.folder, SCHYPE)
    
   
    # make the motifstype file of the ALL group
    motifs_type_all(args.folder, SCHYPE)
    
    # concatenate files
    with open(f"{args.folder}/{SCHYPE}/SCHypeALL/2nodeType.txt", "wb") as Output:
        shutil.copyfileobj(open(f"{args.folder}/2nodemotifs/DD2nodemotifs.txt", "rb"), Output)
        shutil.copyfileobj(open(f"{args.folder}/2nodemotifs/DU2nodemotifs.txt", "rb"), Output)

    # cut files
    with open(f"{args.folder}/{SCHYPE}/SCHypeALL/2nodeType.txt", "r") as f_type:
        with open(f"{args.folder}/{SCHYPE}/SCHypeALL/2node.txt", "w") as f:
            for line in f_type:
                line = line.rstrip()
                line = line.split('\t')
                f.write(f"{line[0]}\t{line[1]}\n")

    # cut files
    with open(f"{args.folder}/{SCHYPE}/SCHypeALL/MotifsType.txt", "r") as f_type:
        with open(f"{args.folder}/{SCHYPE}/SCHypeALL/Motifs.txt", "w") as f:
            for line in f_type:
                line = line.rstrip()
                line = line.split("\t")
                f.write(f"{line[0]}\t{line[1]}\t{line[2]}\n")

    # concatenate files
    with open(f"{args.folder}/{SCHYPE}/SCHypeALL/MotifsALL.txt", "wb") as Output:
        shutil.copyfileobj(open(f"{args.folder}/{SCHYPE}/SCHypeALL/Motifs.txt", "rb"), Output)
        shutil.copyfileobj(open(f"{args.folder}/{SCHYPE}/SCHypeALL/2node.txt", "rb"), Output)
    
    # make the motif files (input files SCHype) for each group
    motif_file(args.folder, Type_hash, SCHYPE)
    
    # collect all types of subgraphs (unique values)
    Types = ['ALL']
    Types += list(set([x for x in Type_hash]))

    # for each group, execute the SCHype
    for type_motif in Types:
        execute(args.folder, SCHYPE, type_motif, args.Path_to_SCHYPE_jar, args.P_value_SCHYPE)
    
    # execute SCHype for each file in the ISMA folder
    # execute_momo(args.folder, SCHYPE, args.Path_to_SCHYPE_jar)

    try:
        os.mkdir(f"{args.folder}/Logs")
    except FileExistsError:
        pass
    with open(f"{args.folder}/Logs/SCHYPE_done.txt", 'w+') as f:
        f.write("SCHYPE done!")
