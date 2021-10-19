import sys
import os

file    = sys.argv[1]
species = sys.argv[2]
ui      = sys.argv[3].split('_') 
di      = sys.argv[4].split('_')

try:
    os.mkdir(f'{species}')
except FileExistsError:
    pass

try:
    os.mkdir(f'{species}/Interactions')
except FileExistsError:
    pass

#ui = ['P', 'H', 'G']                    # Contains all undirected interactions
#di = ['D', 'd', 'R', 'r', 'M', 'm']     # Contains all directed interactions

#print(ui)
#print(di)

with open(file, 'r') as f:
    for line in f:
        line = line.rstrip()
        id1, id2, interaction = line.split('\t')
        if interaction in ui:
            direction = 'u'
            with open(f'{species}/Interactions/{species}_{interaction.upper()}_{direction}.txt', 'a+') as i:
                i.write(f'{id1}\t{id2}\n')
        else:
            direction = 'd'
            interaction_inv = interaction.lower()
            with open(f'{species}/Interactions/{species}_{interaction.upper()}_{direction}.txt', 'a+') as i:
                    #open(f'{species}/Interactions/{species}reverse_{interaction_inv}_{direction}.txt', 'a+') as j:
                i.write(f'{id1}\t{id2}\n')
                #j.write(f'{id2}\t{id1}\n')


try:
    os.mkdir(f"{species}/Logs")
except FileExistsError:
    pass

with open(f"{species}/Logs/split_interactions_done.txt", 'w+') as f:
    f.write("splitting interactions done!")
