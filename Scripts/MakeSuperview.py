"""
Created on 17-okt.-2016
@author: jofoo
"""
import csv
import pandas as pd
import itertools
import os
import sys

# read interactions

Folder = sys.argv[1]

Groups = [G.replace('SCHype', '') for G in os.listdir(f'{Folder}/SCHYPE')]  # UNDEF, MOMO, CIR, COR, COP, COM, FBU, FB2U, FFL, ALL
if not os.path.exists(f'{Folder}/Superview'):
    os.makedirs(f'{Folder}/Superview')

interactions = [f'{Folder}/Interactions/{G}' for G in os.listdir(f'{Folder}/Interactions')]
interaction_dict = {}
for data in interactions:
    if data[-7].isupper():
        interaction_dict[data[-7] + 'int'] = pd.read_csv(data, sep='\t', header=0, names=['A', 'B'])

for SCHypeGroup in Groups:
    NewNodes = {}  # translate modules to one node with gene list as feature
    with open((Folder + '/SCHYPE/SCHype' + SCHypeGroup + '/' + SCHypeGroup + '.nodes.txt'), 'r') as nodefile:
        nodescluster = csv.DictReader(nodefile, delimiter='\t', dialect="excel", fieldnames=['gene', 'cluster'])
        for row in nodescluster:
            if row['cluster'] in NewNodes.keys():
                NewNodes[row['cluster']] = NewNodes[row['cluster']] + '|' + row['gene']
            else:
                NewNodes[row['cluster']] = row['gene']

    # filter cluster nodes between 5 and 50
    deletarray = []
    ClusterSize = {}
    for Cl in NewNodes.keys():
        Genes = NewNodes[Cl].split('|')
        ClusterSize[Cl] = len(Genes)
        if len(Genes) < 5:
            deletarray.append(Cl)
        elif len(Genes) > 50:
            deletarray.append(Cl)

    # delete small/big cluster
    print(SCHypeGroup + " sizefilter: " + str(len(deletarray)))
    for ElementToDelete in deletarray:
        del NewNodes[ElementToDelete]

    # Delete full homolog clusters
    Hcount = {}
    EdgeCount = {}
    with open((Folder + '/SCHYPE/SCHype' + SCHypeGroup + '/ALLEDGES.txt'), 'r') as edgefile:
        edgecluster = csv.DictReader(edgefile, delimiter='\t', dialect="excel",
                                     fieldnames=['cluster', 'gene1', 'interaction', 'gene2'])
        for row in edgecluster:
            if row['cluster'] in Hcount.keys() and row['interaction'] == 'H':
                Hcount[row['cluster']] = Hcount[row['cluster']] + 1
            elif row['interaction'] == 'H':
                Hcount[row['cluster']] = 1

            if row['cluster'] in EdgeCount.keys():
                EdgeCount[row['cluster']] = EdgeCount[row['cluster']] + 1
            else:
                EdgeCount[row['cluster']] = 1

    deletarray = []
    for Cl in Hcount.keys():
        if Hcount[Cl] / EdgeCount[Cl] > 0.9:
            deletarray.append(Cl)

    # delete small/big cluster
    print(SCHypeGroup + " Hfilter: " + str(len(deletarray)))
    for ElementToDelete in deletarray:
        if ElementToDelete in NewNodes.keys():
            del NewNodes[ElementToDelete]

    # Merge clusters which are 50 % overlapping
    def intersect(a, b):
        return list(set(a) & set(b))

    deletarray = []
    for Cl in NewNodes.keys():
        Genes = NewNodes[Cl].split('|')
        for Cl2 in NewNodes.keys():
            if Cl != Cl2:
                Genes2 = NewNodes[Cl2].split('|')
                if len(intersect(Genes, Genes2)) / len(list(set(Genes + Genes2))) > 0.5:

                    # merge small one in bigger one, delete te small one
                    if len(Genes) >= len(Genes2):
                        NewNodes[Cl] = '|'.join(list(set(Genes + Genes2)))
                        if Cl2 not in deletarray:
                            deletarray.append(Cl2)
                    else:
                        NewNodes[Cl2] = '|'.join(list(set(Genes + Genes2)))
                        if Cl not in deletarray:
                            deletarray.append(Cl)

    print(SCHypeGroup + " merged: " + str(len(deletarray)))
    print(deletarray)
    for ElementToDelete in deletarray:
        del NewNodes[ElementToDelete]

    GenesToClusters = {}  # per gene in which cluster they occur
    for Cl in NewNodes.keys():
        Genes = NewNodes[Cl].split('|')
        for G in Genes:
            if G in GenesToClusters.keys():
                GenesToClusters[G] = GenesToClusters[G] + '|' + Cl
            else:
                GenesToClusters[G] = Cl

    # Create superview.txt for each interaction type
    def superview(int_types, int_data):
        int_out = pd.DataFrame()
        for n in range(0, len(int_data) - 1):
            if (int_data.A[n] in GenesToClusters.keys()) and (int_data.B[n] in GenesToClusters.keys()):  # check if interaction is clustered
                clusters1 = GenesToClusters[int_data.A[n]].split('|')
                clusters2 = GenesToClusters[int_data.B[n]].split('|')
                if clusters1 != clusters2:
                    int_out = int_out.append(pd.DataFrame(list(itertools.product(clusters1, clusters2))))
        if int_out.empty:
            return int_out
        int_out.values.sort()  # filter so the smallest clusternumber is fist colum
        int_out = int_out[int_out.apply(lambda x: min(x) != max(x), 1)]  # filter out within modules
        int_count_matrix = int_out.groupby([0, 1]).size().reset_index(name='count')

        int_count_matrix[0] = str(SCHypeGroup + '_') + int_count_matrix[0].astype(str)
        int_count_matrix[1] = str(SCHypeGroup + '_') + int_count_matrix[1].astype(str)
        int_count_matrix['C'] = int_types[0]
        int_count_matrix = int_count_matrix[[0, 'C', 1, 'count']]
        return int_count_matrix


    for int_type in interaction_dict.keys():
        int_countMatrix = superview(int_type, interaction_dict[int_type])
        int_countMatrix.to_csv((Folder + '/Superview/SuperView_' + int_type[0] + '_' + SCHypeGroup + '.txt'), sep='\t',
                               header=True, index=False, mode='w')
        print(SCHypeGroup, int_type, 'Superview Created')

    # TFs vs modules
    if 'Dint' in interaction_dict:
        TFout = pd.DataFrame()
        TFlist = list(set(interaction_dict['Dint'].A))
        for TF in range(0, len(TFlist) - 1):
            TFtargets = list(interaction_dict['Dint'].loc[interaction_dict['Dint']['A'] == TFlist[TF], 'B'])
            for Cl in NewNodes.keys():
                Genes = NewNodes[Cl].split('|')
                if TFlist[TF] not in Genes:
                    if len(intersect(Genes, TFtargets)) > 0:
                        TFout = TFout.append(
                            pd.DataFrame({'A': TFlist[TF], 'B': Cl, 'count': len(intersect(Genes, TFtargets))},
                                         index=['a']))

        TFcountMatrix = TFout
        if not TFcountMatrix.empty:
            TFcountMatrix.reset_index(inplace=True)
            TFcountMatrix['B'] = str(SCHypeGroup + '_') + TFcountMatrix['B'].astype(str)
            TFcountMatrix['C'] = "D"
            TFcountMatrix = TFcountMatrix[['A', 'C', 'B', 'count']]
        TFcountMatrix.to_csv((Folder + '/Superview/SuperView_TF_' + SCHypeGroup + '.txt'), sep='\t', header=True,
                             index=False, mode='w')
        print(SCHypeGroup, 'TF Superview Created')

    # miRNA vs modules
    Mout = pd.DataFrame()
    Mlist = list(set(interaction_dict['Mint'].A))
    for M in range(0, len(Mlist) - 1):
        Mtargets = list(interaction_dict['Mint'].loc[interaction_dict['Mint']['A'] == Mlist[M], 'B'])
        for Cl in NewNodes.keys():
            Genes = NewNodes[Cl].split('|')
            if Mlist[M] not in Genes:
                if len(intersect(Genes, Mtargets)) > 0:
                    Mout = Mout.append(
                        pd.DataFrame({'A': Mlist[M], 'B': Cl, 'count': len(intersect(Genes, Mtargets))}, index=['a']))

    McountMatrix = Mout
    if not McountMatrix.empty:
        McountMatrix.reset_index(inplace=True)
        McountMatrix['B'] = str(SCHypeGroup + '_') + McountMatrix['B'].astype(str)
        McountMatrix['C'] = "M"
        McountMatrix = McountMatrix[['A', 'C', 'B', 'count']]
    McountMatrix.to_csv((Folder + '/Superview/SuperView_M_' + SCHypeGroup + '.txt'), sep='\t', header=True, index=False,
                        mode='w')
    print(SCHypeGroup, 'M Superview Created')

    # Get all links between the new modules
    Iout = pd.DataFrame()
    for Mod in NewNodes.keys():
        Genes = NewNodes[Mod].split('|')
        for G in Genes:
            Clusters = GenesToClusters[G].split('|')
            # I = internal = same node in 2 modules
            for C in Clusters:
                if C != Mod:
                    tmp = [Mod, C]
                    tmp.sort()
                    df2 = pd.DataFrame([tmp])
                    Iout = Iout.append(df2)
    if Iout.empty:
        IcountMatrix = Iout
    else:
        Iout.values.sort()  # filter so the smallest clusternumber is fist colum
        Iout = Iout[Iout.apply(lambda x: min(x) != max(x), 1)]  # filter out within modules
        IcountMatrix = Iout.groupby([0, 1]).size().reset_index(name='count')
        IcountMatrix[0] = str(SCHypeGroup + '_') + IcountMatrix[0].astype(str)
        IcountMatrix[1] = str(SCHypeGroup + '_') + IcountMatrix[1].astype(str)
        IcountMatrix['C'] = "I"
        IcountMatrix = IcountMatrix[[0, 'C', 1, 'count']]
    IcountMatrix.to_csv((Folder + '/Superview/SuperView_I_' + SCHypeGroup + '.txt'), sep='\t', header=True, index=False,
                        mode='w')
    print(SCHypeGroup, 'I Created')


try:
    os.mkdir("Logs")
except FileExistsError:
    pass
with open("Logs/Super_View_done.txt", 'w+') as f:
    f.write("Super view done!")
