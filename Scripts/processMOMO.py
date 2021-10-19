import sys
import os

RUN = sys.argv[1]       # location of the MOMO group output of SCHYPE

sub_graph_files = []
with open(f'{RUN}/Motiffiles.txt', 'r') as f:
    for row in f:
        row = row.rstrip().replace('.txt', '').replace('-', '')
        sub_graph_files.append(row)

number = 0
extra = 0
node_file = []
for file in sub_graph_files:
    with open(f'{RUN}/{file}.nodes.txt', 'r') as f:
        for line in f:
            line = line.rstrip()
            Split = line.split('\t')
            if Split[1] == number:
                nrtoPrint = number + extra
                node_file.append(f"{Split[0]}\t{nrtoPrint}")
            else:
                number += 1
                nrtoPrint = number + extra
                node_file.append(f"{Split[0]}\t{nrtoPrint}")
        extra += number + 1
        number = 0


number = 0
extra = 0
edge_file = []
sub_graph_type_file = []
for file in sub_graph_files:
    with open(f"{RUN}/{file}.edges.txt", 'r') as f:
        for line in f:
            line = line.rstrip()
            Split = line.split('\t')
            types = ''.join(i for i in file if not i.isdigit())  # remove digits in file
            if Split[0] == number:
                nrtoPrint = number + extra
                edge_file.append(f"{nrtoPrint}\t{Split[1]}\t{Split[2]}\t{Split[3]}")
                sub_graph_type_file.append(f"{Split[1]}\t{Split[2]}\t{Split[3]}\t{types}")

            else:
                number += 1
                nrtoPrint = number + extra
                edge_file.append(f"{nrtoPrint}\t{Split[1]}\t{Split[2]}\t{Split[3]}")
                sub_graph_type_file.append(f"{Split[1]}\t{Split[2]}\t{Split[3]}\t{types}")
        extra += number + 1
        number = 0


with open(f'{RUN}/MOMO.nodes.txt', 'w+') as f:
    for node in node_file:
        f.write(node + '\n')

with open(f'{RUN}/MOMO.edges.txt', 'w+') as f:
    for edge in edge_file:
        f.write(edge + '\n')

with open(f'{RUN}/MotifsType.txt', 'w+') as f:
    for motiftype in sub_graph_type_file:
        f.write(motiftype + '\n')

try:
    os.mkdir("Logs")
except FileExistsError:
    pass
with open("Logs/process_MOMO_done.txt", 'w+') as f:
    f.write("MOMO processing done!")
