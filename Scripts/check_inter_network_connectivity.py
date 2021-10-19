import argparse
import seaborn
import pandas as pd

def arg_parse():
    """
    Parses the arguments given in the command line
    """
    parser = argparse.ArgumentParser(description='Check out the interaction file and make a statistic on how many edges are shared between network types.')
    parser.add_argument('-interactions', help='Interaction_file.')
    parser.add_argument('-output', help='Path and name to save output')

    return parser.parse_args()


def loadInteractions(interactions):
    """
    Load the interaction file
    """
    new_sets = {}
    new_list = {}

    infile = open(interactions, 'r')
    lines  = infile.read().split('\n')

    for line in lines:
        tmp = line.split('\t')
        if len(line)>1:
            if tmp[2] not in new_sets:
                new_sets[tmp[2]] = {tmp[0]:'', tmp[1]:''}
                new_list[tmp[2]] = []
            else:
                new_sets[tmp[2]][tmp[0]] = ''
                new_sets[tmp[2]][tmp[1]] = ''
            new_list[tmp[2]].append([tmp[0], tmp[1]])

    return new_sets, new_list










if __name__ == '__main__':
    """
    Check out the interaction file and make a statistic on how many edges are shared between
    network types
    """

    #parse argumnets
    args = arg_parse()

    # load interactions
    interactions_set, interactions_list = loadInteractions(args.interactions)

    undirected = ["P", "H", "G"]

    df_store = []

        

    #print(pd.DataFrame([[1, 2, 3]], columns=(["Group", "ID", "Values"])))
    #exit()
 
    # for all network edes
    for k1 in interactions_list.keys():
       # against all network edges
        store = [k1]
        for k2 in interactions_list.keys(): 
            counter = 0
            for n in interactions_list[k1]:
                if n[0] in interactions_set[k2] or n[1] in interactions_set[k2]:
                    for k in interactions_list[k2]:
                        if n[0] == k[0]:
                            counter += 1
                        if n[1] == k[0]:
                            counter += 1
                        if n[1] == k[0]:
                            counter += 1
                        if n[1] == k[1]:
                            counter += 1
            #print(k1, k2, len(interactions_list[k1]), len(interactions_list[k2]), counter)
            store.append(counter)
            df_store.append(pd.DataFrame([[k1, k2, counter]], columns=(["Group", "ID", "Values"])))


    df = df_store[0]
    for i in df_store[1:]:
        df.append(i, ignore_index=True)
    print(df)
    fg = seaborn.factorplot(x='ID', y='Values', hue='Group', kind='bar', data=df)
    print(df)


