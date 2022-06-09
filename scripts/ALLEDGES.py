import os
import numpy as np
import pandas as pd
import itertools
import sys

Folder = sys.argv[1]
Groups = [G.replace('SCHype', '') for G in os.listdir(Folder + '/SCHYPE')]

for SCHypeGroup in Groups:
    # Create reference dictionary
    with open((Folder + '/SCHYPE/SCHype' + SCHypeGroup + '/MotifsType.txt'), 'r') as SubgraphsType:
        interaction_dict = {}
        types = ''
        tuple_list = []

        for row in SubgraphsType:
            row = row.rstrip().split('\t')
            if row[3] == types:
                tuple_list.append(tuple(row[0:3]))
            else:
                if types != '':
                    interaction_dict[types] = tuple_list
                types = row[3]
                tuple_list = [tuple(row[0:3])]
        interaction_dict[types] = tuple_list

    # Match edges with their interaction type
    with open((Folder + '/SCHYPE/SCHype' + SCHypeGroup + '/' + SCHypeGroup + '.edges.txt'), 'r') as edgefile:
        cluster_num = '0'
        edge_searched_list = []
        ALLEDGES = pd.DataFrame()
        for row in edgefile:
            row = row.rstrip().split('\t')
            if cluster_num != row[0]:  # Refresh edge_list for each cluster
                cluster_num = row[0]
                edge_searched_list = []

            edge_tuple = tuple(row[1:])
            for interaction_type in interaction_dict.keys():
                if edge_tuple in interaction_dict[interaction_type]:
                    comb = list(itertools.combinations(row[1:], 2))
                    for edge in range(3):
                        if comb[edge] not in edge_searched_list:
                            edge_searched_list.append(comb[edge])
                            ALLEDGES = ALLEDGES.append(pd.DataFrame(
                                np.array([[row[0], comb[edge][0], interaction_type[edge], comb[edge][1]]]),
                                index=['a']))
                    break
    # Write ALLEDGES.txt
    ALLEDGES.to_csv((Folder + '/SCHYPE/SCHype' + SCHypeGroup + '/' + 'ALLEDGES.txt'), sep='\t', header=False,
                    index=False, mode='w')
    print(SCHypeGroup, 'Finished')

#print(sys.argv)
try:
    os.mkdir(sys.argv[1]+"/Logs")
except FileExistsError:
    pass
with open(sys.argv[1]+"/Logs/ALLEDGES_done.txt", 'w+') as f:
    f.write("alledges.txt done!")
