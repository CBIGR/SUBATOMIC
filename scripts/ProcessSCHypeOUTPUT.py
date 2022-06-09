import argparse
import os
import re
import shutil


def merge_subgraph(subgraphs):
    """
    Merges the subgraph types into general type.
        E.g.: Merged_subgraph = [D,D-P,D-M]
    """

    # joins together the first, second and third interaction of the subgraph with a '-'
    merged_subgraph = ['-'.join([m[0] for m in subgraphs]),
                       '-'.join([m[1] for m in subgraphs]),
                       '-'.join([m[2] for m in subgraphs])]

    return merged_subgraph


def motif_to_edges(edges, subgraph, di, ui):
    """
    Reverts the edges of the subgraphs with lower case.
    Removes any duplicates in undirected edges.
    Returns the cluster and a unique set of the modified edges.
    """
    edges_out, edges_out_n = [], []

    # Motif with lower case revert edge
    edges_in = edges.split('\t')
    cluster = edges_in[0]

    # First edge
    if '-' not in subgraph[0]:
        if subgraph[0].isupper():               # Test if uppercase
            edges_out.append(f'{edges_in[1]}\t{subgraph[0]}\t{edges_in[2]}')
        else:                                   # Lower case revert edges
            edges_out.append(f'{edges_in[2]}\t{subgraph[0].upper()}\t{edges_in[1]}')
    else:
        types = subgraph[0].split('-')
        for letter in types:
            if letter.isupper():                # Test if uppercase
                edges_out.append(f'{edges_in[1]}\t{letter}\t{edges_in[2]}')
            else:                               # Lower case revert edges
                edges_out.append(f'{edges_in[2]}\t{letter.upper()}\t{edges_in[1]}')

    # Second edge
    if '-' not in subgraph[1]:
        if subgraph[1].isupper():               # Test if uppercase
            edges_out.append(f'{edges_in[1]}\t{subgraph[1]}\t{edges_in[3]}')
        else:                                   # Lower case revert edges
            edges_out.append(f'{edges_in[3]}\t{subgraph[1].upper()}\t{edges_in[1]}')
    else:
        types = subgraph[1].split('-')
        for letter in types:
            if letter.isupper():                # Test if uppercase
                edges_out.append(f'{edges_in[1]}\t{letter}\t{edges_in[3]}')
            else:                               # Lower case revert edges
                edges_out.append(f'{edges_in[3]}\t{letter.upper()}\t{edges_in[1]}')

    # Third edge
    if '-' not in subgraph[2]:
        if subgraph[2].isupper():               # Test if uppercase
            edges_out.append(f'{edges_in[2]}\t{subgraph[2]}\t{edges_in[3]}')
        else:                                   # Lower case revert edges
            edges_out.append(f'{edges_in[3]}\t{subgraph[2].upper()}\t{edges_in[2]}')
    else:
        types = subgraph[2].split('-')
        for letter in types:
            if letter.isupper():                # Test if uppercase
                edges_out.append(f'{edges_in[2]}\t{letter}\t{edges_in[3]}')
            else:                               # Lower case revert edges
                edges_out.append(f'{edges_in[3]}\t{letter.upper()}\t{edges_in[2]}')


    return cluster, list(set(edges_out))


def functional_description(file):
    """
    Splits up the functional description file.
    Functional description file contains on each line:
        gene id
        public name
        molecular name
        concise description
        TF type
    All separated by a tab, in that specific order.
    For unknown values, the term "not known" should be written.
    Lines that should be ignored should start with '#'.
    Returns:
        A dictionary containing the accession number as keys and a name for the gene, if specified, as values;
        A dictionary containing the gene ID as values and the public name as keys;
        A dictionary containing the locus ID with the molecular name as keys and the gene id's as values;
        A list with all the functional descriptions.
    """
    accession_nr, id_gene, id_locus, biotype = {}, {}, {}, {}
    func_names = []

    with open(file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line[0] != '#':          # ignore lines starting with '#'

                # split up the line in the stated parts
                gene_id, public_name, molecular_name, concise_description, biotype_single = line.split('\t')[:5]

                # add the different parts to the corresponding dictionary and to the list with functional descriptions
                if public_name.lower() != 'not_known':
                    accession_nr[gene_id] = public_name
                    func_names.append(f'{gene_id}\t{public_name}\t{concise_description}')
                elif molecular_name.lower() != 'not_known':
                    accession_nr[gene_id] = molecular_name
                    func_names.append(f'{gene_id}\t{molecular_name}\t{concise_description}')
                else:
                    accession_nr[gene_id] = gene_id
                    func_names.append(f'{gene_id}\t{gene_id}\t{concise_description}')

                if public_name != '':
                    id_gene[public_name] = gene_id
                if molecular_name != '':
                    id_locus[molecular_name] = gene_id

                biotype[gene_id] =  biotype_single
                #print(biotype_single)

                

    return accession_nr, id_gene, id_locus, func_names, biotype


def transcription_factors(tf_file):
    """
    Splits up the transcription factor file.
    One line in the transcription file contains:
        The name of the transcription factor;
        The type of the transcription factor.
    Separated by a tab, in that specific order.
    Unknown values can be left blank.
    Returns dictionaries of the gene type and the TF type.
    """

    gene_type, tf_type = {}, {}
    with open(tf_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            tf = line.split('\t')

            # check if the TF type is known
            if len(tf) == 2:
                type_tf = line.split('\t')[1]
            else:
                type_tf = 'not_known'

            # add the TF to the gene type dictionary
            gene_type[tf[0]] = 'TF'

            # add the TF type to the TF type dictionary
            # not_known will be added if the type is unknown
            tf_type[tf[0]] = type_tf

    return gene_type, tf_type


def mirna(mirna_file, gene_type):
    """
    Splits up the miRNA file.
    One line in the miRNA file contains:
        gene id of the gene which is a miRNA.
    Returns dictionary of the gene type.
    """

    with open(mirna_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            gene_type[line] = 'miRNA'       # adds the miRNA genes to the gene type dictionary

    return gene_type


def groups(module_data):
    """
    Returns a list of all the groups present in the given folder of the module data, sorted alphabetically
    """

    return sorted([x for x in os.walk(module_data)][0][1])


def gene_ontology(go_file, accession_nr, id_gene):
    """
    Reads in the properties of the genes with GO.
    Returns a dictionary containing the properties of the genes with GO.
    If no GO file is provided, None is returned.
    """
    go = {}
    if go_file is None:
        go = None
    else:
        with open(go_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                go[line] = "YES"
                go[id_gene[line]] = "YES"
                go[accession_nr[line]] = "YES"

    return go


def merge_node_attributes(func_names, gene_type, tf_type, ontology):
    """
    Merges the node attributes.
    Returns a list containing the node attributes.
        Node attributes:
            accession number of the gene;
            name of the gene, if given in the functional description file;
            the type of the gene, if not specified, the gene type is just gene;
            if the gene is a TF, the type of TF, if this was given in the TF file;
            the GO of the gene, if given in the GO file.
        All these are written on one line, separated by a tab.
    """
    node_attributes = []
    for func in func_names:
        func = func.rstrip()
        gene_id, gene_name, gene_description = func.split('\t')

        if gene_id in gene_type:
            attributes = f'{gene_id}\t{gene_name}\t{gene_type[gene_id]}'
        else:
            attributes = f'{gene_id}\t{gene_name}\tGene'

        if gene_id in tf_type:
            attributes += f'\t{tf_type[gene_id]}'
        else:
            attributes += '\t-'

        if ontology:
            if gene_id in ontology:
                attributes += f'\t{ontology[gene_id]}'
        else:
            attributes += '\t-'

        node_attributes.append(attributes)

    return node_attributes


def two_node_hash(subgraph_type):
    """
    Makes a dictionary:
        Keys: the two nodes, separated by a tab;
        Values: the subgraph
    The input is the 2nodeType.txt file created by the prepSCHype.py script:
        On each line: node one \t node two \t subgraph
    """
    type_hash = {}
    with open(subgraph_type, 'r') as f:
        for node in f:
            node = node.rstrip()
            node_one, node_two, subgraph = node.split('\t')

            if f'{node_one}\t{node_two}' not in type_hash:
                type_hash[f'{node_one}\t{node_two}'] = []
            type_hash[f'{node_one}\t{node_two}'].append(subgraph)

    return type_hash


def three_node_hash(subgraph_type):
    """
    Makes a dictionary:
        Keys: the three nodes, separated by a tab;
        Values: the subgraph
    The input is the MotifsType.txt file created by the prepSCHype.py script:
        On each line: node one \t node two \t node three \t subgraph
    """
    type_hash = {}
    with open(subgraph_type, 'r') as f:
        for line in f:
            line = line.rstrip()
            node_one, node_two, node_three, motif = line.split('\t')

            if f'{node_one}\t{node_two}\t{node_three}' not in type_hash:
                type_hash[f'{node_one}\t{node_two}\t{node_three}'] = []
            type_hash[f'{node_one}\t{node_two}\t{node_three}'].append(motif)

    return type_hash


def two_node_subgraph(edges, node_2_type_hash, sif, nnf, eda, clusters, group):
    """
    The edges are coming from the group.edges.txt file, created by SCHype.
    Adds the edges, given by the node_2_type_hash (made earlier) to the sif, nnf, eda and clusters, which are used
    later to create the corresponding files. Does this for the two node subgraphs.
    """

    cluster, node_one, node_two = edges
    for motif in node_2_type_hash[f'{node_one}\t{node_two}']:

        # make the edges in both directions
        edge_one = f'{node_one}\t{motif[0]}\t{node_two}'
        edge_two = f'{node_two}\t{motif[1]}\t{node_one}'

        # add the first edge to the sif, nnf, clusters and eda
        sif.append(edge_one)
        nnf.append(f'{group[6:]}_{cluster}\t{edge_one}')
        clusters.append(cluster)
        edge_one = edge_one.replace('\t', ' (', 1)
        edge_one = edge_one.replace('\t', ') ', 1)
        eda.append(f'{edge_one}\t{group[6:]}_{cluster}')

        # add the second edge to the sif, nnf, clusters and eda
        sif.append(edge_two)
        nnf.append(f'{group[6:]}_{cluster}\t{edge_two}')
        clusters.append(cluster)
        edge_two = edge_two.replace('\t', ' (', 1)
        edge_two = edge_two.replace('\t', ') ', 1)
        eda.append(f'{edge_two}\t{group[6:]}_{cluster}')

    return sif, nnf, eda, clusters


def three_node_subgraph(edges, node_3_type_hash, line, sif, nnf, eda, clusters, cluster_edges, group, di, ui):
    """
    The edges are coming from the group.edges.txt file, created by SCHype.
    Adds the edges, given by the node_3_type_hash (made earlier) to the sif, nnf, eda, clusters and cluster edges,
    which are used later to create the corresponding files. Does this for the three node subgraphs.
    """

    # subgraphs coming from the node_3_type_hash, made earlier
    subgraphs = node_3_type_hash[f'{edges[1]}\t{edges[2]}\t{edges[3]}']
    
    # Merge the subgraph types into general type, e.g. Merged motif [D,D-P,D-M]
    merged_motif = merge_subgraph(subgraphs)
   
    # Sub input line return list with clusters and list with edges
    cluster, edges_out = motif_to_edges(line, merged_motif, di, ui)
        
    for edge in edges_out:
        # Write the edges to appropriate lists in right format
        sif.append(edge)
        nnf.append(f'{group[6:]}_{cluster}\t{edge}')
        clusters.append(cluster)
        edge = edge.replace('\t', ' (', 1)
        edge = edge.replace('\t', ') ', 1)
        eda.append(f'{edge}\t{group[6:]}_{cluster}')

        edge = re.sub(r'\(.*\)', '', edge)
        if f'{group[6:]}_{cluster}' not in cluster_edges:
            cluster_edges[f'{group[6:]}_{cluster}'] = ''
        cluster_edges[f'{group[6:]}_{cluster}'] += f'{edge}\n'

    return sif, nnf, eda, clusters, cluster_edges


def clustering(edges_schype, group, two_nodes_hash, three_nodes_hash, di, ui):
    """
    Creates the lists sif, nnf, eda, cluster and the dictionary cluster edges, which are later used to make the
    corresponding files. Does this for the two node subgraphs and the three node subgraphs.
    """
    cluster_edges = {}
    sif, nnf, eda, cluster = [], [], [], []

    with open(edges_schype, 'r') as f:
        for line in f:
            line = line.rstrip()
            edges = line.split('\t')

            # 2-node subgraph
            if group == 'SCHypeALL' and len(edges) == 3:
                sif, nnf, eda, cluster = two_node_subgraph(edges, two_nodes_hash, sif, nnf, eda, cluster, group)

            # 3-node subgraph
            else:
                sif, nnf, eda, cluster, cluster_edges = three_node_subgraph(edges, three_nodes_hash, line, sif, nnf,
                                                                            eda, cluster, cluster_edges, group, di, ui)

    return sif, nnf, cluster, eda, cluster_edges


def node_file(file, group):
    """
    Splits up the group.nodes.txt file, created by SCHype.
    Returns a dictionary:
        Keys: group_cluster;
        Values: gene.
    """

    cluster_genes = {}

    with open(file, 'r') as f:
        for line in f:
            line = line.rstrip()
            gene, cluster = line.split('\t')

            if f'{group}_{cluster}' not in cluster_genes:
                cluster_genes[f'{group}_{cluster}'] = []
            cluster_genes[f'{group}_{cluster}'].append(f'{gene}\t')

    return cluster_genes


def cluster_size_attributes(cluster_genes, tf_type, gene_go):
    """
    Calculates the number of genes with GO in a cluster.
    Calculates the number of transcription factors in a cluster.
    Returns this as a dictionary (cluster_attributes):
        Keys: Cluster;
        Values: List: [TF count, GO count, cluster size].
    Returns a dictionary (cluster_size):
        Keys: Cluster;
        Values: cluster size.
    """

    cluster_size, cluster_attributes = {}, {}
    for cluster in cluster_genes:
        genes = cluster_genes[cluster]
        cluster_size[cluster] = len(genes)
        go_count = 0
        tf_count = 0
        for gene in genes:
            if gene.replace('\t', "") in tf_type:
                tf_count += 1
            if gene_go and gene.replace('\t', "") in gene_go:
                go_count += 1
        
        if cluster not in cluster_attributes:
            cluster_attributes[cluster] = [0, 0, 0]
        cluster_attributes[cluster][0] += tf_count
        cluster_attributes[cluster][1] += go_count
        cluster_attributes[cluster][2] += len(genes)

    return cluster_size, cluster_attributes


def module_viewer(base_folder, run, group, cluster_genes):
    """
    Writes the module viewer file:
        On each line:
            Cluster \t the genes from that cluster separated by |.
    """

    for cluster in cluster_genes:
        genes = cluster_genes[cluster]
        if len(genes) > 5:
            with open(f'{base_folder}/{run}/ModuleViewer/data/{group[6:]}.txt', 'a') as f:
                f.write(f'{cluster[6:]}\t{"|".join(genes)}\n')


def sif_file(folder, group, sif):
    """
    Makes SIF file for visualization in Cytoscape.
    Based on the sif list made earlier.
    On each line:
        GENE1 \t Edge type \t GENE2
    """
    with open(f'{folder}/{group}/edges.sif', 'w') as f:
        for s in sif:
            f.write(s + '\n')


def nnf_file(folder, group, clusters, nnf, hom):
    """
    Makes the nnf file for visualization in Cytoscape.
    Made based on the nnf list made earlier.
    On each line:
        cluster \t GENE1 \t EdgeType \t GENE2
    Returns:
        Edge value: number of edges in a cluster;
        Hom value: the number of homologous interactions in a cluster.
    """
    edge_value, hom_value = {}, {}

    # check wich network types contain home or miRNA information and add interaction
    # letters to variables
    hm_networks = []
    mi_networks = []
    if hom != None:
        hm_networks = hom.split('_')

    lists = {}
    with open(f'{folder}/{group}/edges.nnf', 'w') as f:
        for cluster in sorted(list(set(clusters))):
            f.write(f'Network\t{group[6:]}_{cluster}\n')

        for n in sorted(list(set(nnf))):
            f.write(f'{n}\n')               # cluster \t gene1 \t edge type \t gene2
            cluster, gene_one, edge_type, gene_two = n.split('\t')
            lists[cluster] = ''

            # count the number of edges in a cluster
            if cluster not in edge_value:
                edge_value[cluster] = 0
            edge_value[cluster] += 1

            # count the number of homologous interactions in a cluster
            if edge_type in hm_networks:
                if cluster not in hom_value:
                    hom_value[cluster] = 0
                hom_value[cluster] += 1

    return edge_value, hom_value


def write_attributes(folder, group, attributes, edge_values, hom_values, attribute_type):
    """
    Makes and (over)writes the attribute files
    "Attributes" contains the specific attribute dictionary
    "Attribute_type" can be:
        Cluster: .cla file Cluster \t Number of edges in cluster \t fraction of homologous interactions in cluster
        Edge: .eda file GENE1 \t (EdgeType) \t GENE2 \t cluster
        Node: .noa file GENE vb: extra names
    """

    # open the file complementary to the attribute type
    if attribute_type.lower() == 'cluster':
        file = f'{folder}/{group}/cluster.cla'
    elif attribute_type.lower() == 'edge':
        file = f'{folder}/{group}/cluster.eda'
    else:
        file = f'{folder}/{group}/NodeAttributes.noa'

    # write the attributes to the opened file
    # attributes were compiled earlier
    with open(file, 'w') as f:
        if attribute_type.lower() == 'cluster':
            f.write('Cluster\tNrEdges\tPerHOM\n')

            # write the edge value and the hom value if given
            for edge_value in edge_values:
                if edge_value in hom_values:
                    f.write(f'{edge_value}\t{edge_values[edge_value]}\t{round(hom_values[edge_value] / edge_values[edge_value],2)}\n')
                else:
                    f.write(f'{edge_value}\t{edge_values[edge_value]}\t0\n')

        elif attribute_type.lower() == 'edge':
            f.write('Interaction\tCluster\n')
            for attribute in attributes:
                f.write(f'{attribute}\n')

        elif attribute_type.lower() == 'node':
            f.write('Accession_Number\tFuncName\tType\tTF_type\tHas_GO\n')
            for attribute in attributes:
                f.write(f'{attribute}\n')


def data_file(folder, group, clusters, cluster_attributes):
    """
    Makes and (over)writes a datafile per group
    One line in the file contains:
        Cluster
        Comments
        Number of TF in that cluster
        Number of GO genes in that cluster
        nrInCluster
    Separated by a tab and in that specific order
    """
    with open(f'{folder}/{group}/Data_{group[6:]}.txt', 'w') as f:
        f.write('Module\tComments\tnrTF\tnrGO\tnrInCluster\n')
        for a in sorted(list(set(clusters))):
            tmp = f'{group}_{a}'
            if tmp in cluster_attributes:
                tf_count, go_count, nr_of_genes = cluster_attributes[tmp]
                f.write(f'{tmp[6:]}\t-\t{tf_count}\t{go_count}\t{nr_of_genes}\n')


def clustering_coefficient(folder, group, cluster_edges, cluster_coefficient_r, call_r):
    """
    Makes a file with a list of all files containing a cluster.
    Uses that file to calculate the clustering coefficient of that specific cluster.
        Via R file called ClusterCoefficient.R.
    Writes those coefficients to the ClusteringCoefficient.txt file.
    """

    make_dir(f'{folder}/{group}/CC')  # Temporary folder for calculating clustering coefficient

    for cluster_edge in sorted(cluster_edges.keys()):
        with open(f'{folder}/{group}/CC/{cluster_edge}.txt', 'a') as f1:
            with open(f'{folder}/{group}/CC/RCMDlist_{group}.txt', 'a') as f2:
                f1.write(f'{cluster_edges[cluster_edge]}\n')
                f2.write(f'{cluster_edge}.txt\n')

    command = call_r+f' ' \
              f'--vanilla ' \
              f'-f {cluster_coefficient_r} ' \
              f'--args {folder}/{group}/ CC/RCMDlist_{group}.txt'

    os.system(command)

    shutil.rmtree(f"{folder}/{group}/CC")


def arg_parse():
    """
    Parses the command line arguments
    """
    parser = argparse.ArgumentParser(description='Processes the SCHype output')
    parser.add_argument('Module_Data', help='Relative or absolute path to the folder containing the module data')
    parser.add_argument('Output_folder',
                        help='Relative or absolute path to the folder where your output files need to be put')
    parser.add_argument('P_value_SCHYPE', default=1, type=int, help='Value calculates ratio between edges and nodes')
    parser.add_argument('Functional_Description',
                        help='Relative or absolute path to the file containing the functional description')
    parser.add_argument('Transcription_Factors',
                        help='Relative or absolute path to the file containing the transcription factor genes')
    parser.add_argument('R_Clustering_Coefficient',
                        help='Relative or absolute path to the R script calculating the clustering coefficient')
    parser.add_argument('-miRNA', default=None, help='Relative or absolute path to the file containing the miRNA genes')
    parser.add_argument('-GO', default=None, help='Relative or absolute path to the file containing the genes with GO')
    parser.add_argument('-m', default=None, help='String of network letters containing miRNAs')
    parser.add_argument('-hom', default=None, help='String of network letters containing homologs')
    parser.add_argument('-R', default='R', help='How to call R')
    parser.add_argument('-di', default=None, help='Gives network letter from interactions that are directed')
    parser.add_argument('-ui', default=None, help='Gives network letter from interactions that are undirected')


    return parser.parse_args()


def make_dir(name):
    """
    Makes new directories.
    If the directory already exists, it escapes that error.
    """

    try:
        os.mkdir(name)
    except FileExistsError:
        pass

def alterantiveGeneType(biotype):
    """
    Create Gene_Type and TF type from functional annotation file
    """
    
    Gene_Type = {}
    TF_type   = {}

    for i in biotype:
        if biotype[i] != 'Gene':
           Gene_Type[i] = biotype[i]
           TF_type[i] = 'not_known'
    return Gene_Type, TF_type
            
        


if __name__ == '__main__':

    args = arg_parse()

    Base_Folder = args.Output_folder  # Current directory
    SCHYPE = 'SCHYPE'

    #if args.P_value_SCHYPE != 1:
    #    SCHYPE += str(args.P_value_SCHYPE)

    RUN = f'MotifClusters'

    make_dir(f'{Base_Folder}/{RUN}')
    make_dir(f'{Base_Folder}/{RUN}/ModuleViewer')
    make_dir(f'{Base_Folder}/{RUN}/ModuleViewer/data')
    make_dir(f'{Base_Folder}/{RUN}/{SCHYPE}')

    Folder = f'{Base_Folder}/{RUN}/{SCHYPE}'

    Accession_NR, ID_gene, ID_locus, Func_names, biotype = functional_description(args.Functional_Description)  # ID file
    Gene_Type, TF_type = transcription_factors(args.Transcription_Factors)                             # TF file

    #print(Gene_Type, TF_type)
    # create alternative TF type that reads from biotype
    Gene_Type, TF_type = alterantiveGeneType(biotype)
    #print('novel')
    #print(Gene_Type, TF_type)
    #print(len(Gene_Type), len(TF_type))    
    #exit()
    


    if args.miRNA:
        Gene_Type = mirna(args.miRNA, Gene_Type)                                                       # miRNA file

    node_atr = merge_node_attributes(Func_names, Gene_Type, TF_type, gene_ontology(args.GO, Accession_NR, ID_gene))
    



    for Group in groups(args.Module_Data):#["SCHypeCOP"]:#groups(args.Module_Data):
        make_dir(f'{Folder}/{Group}')

        if Group == 'SCHypeALL':
            Two_Node_Hash = two_node_hash(f'{args.Module_Data}/{Group}/2nodeType.txt')
        else:
            Two_Node_Hash = None
        Three_Node_Hash = three_node_hash(f'{args.Module_Data}/{Group}/MotifsType.txt')




        SIF, NNF, Cluster, EDA, Cluster_Edges = clustering(f'{args.Module_Data}/{Group}/{Group[6:]}.edges.txt',
                                                           Group, Two_Node_Hash, Three_Node_Hash, args.di.split('_'), args.ui.split('_'))


        Cluster_Genes = node_file(f'{args.Module_Data}/{Group}/{Group[6:]}.nodes.txt', Group)

        # Could be enhanced by also counting miRNA (?)
        Cluster_Size, Cluster_Attributes = cluster_size_attributes(Cluster_Genes, TF_type, args.GO)


        module_viewer(Base_Folder, RUN, Group, Cluster_Genes)

        sif_file(Folder, Group, SIF)
        EDGE_value, HOM_value = nnf_file(Folder, Group, Cluster, NNF, args.hom)
        
        write_attributes(Folder, Group, Cluster_Attributes, EDGE_value, HOM_value, "cluster")  # cluster
        write_attributes(Folder, Group, EDA, EDGE_value, HOM_value, "edge")                    # edge
        write_attributes(Folder, Group, node_atr, EDGE_value, HOM_value, "node")               # node

        data_file(Folder, Group, Cluster, Cluster_Attributes)
        clustering_coefficient(Folder, Group, Cluster_Edges, args.R_Clustering_Coefficient, args.R)

        print(f'{Group} processed')

    #exit()

    try:
        os.mkdir(f"{Base_Folder}/Logs")
    except FileExistsError:
        pass
    with open(f"{Base_Folder}/Logs/process_SCHYPE_done.txt", 'w') as f:
        f.write("SCHYPE processed!")
