from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Example sequence alignment
alignment_data = {
    'Seq1': 'AGCCG',
    'Seq2': 'AGACG',
    'Seq3': 'ATAAG',
    'Seq4': 'ATGAG',
    'Seq5': 'ACGAG',
    'Seq6': 'ACTAG'
}

# Cluster assignments for each sequence
cluster_assignments = {
    'Seq1': 3,
    'Seq2': 3,
    'Seq3': 2,
    'Seq4': 2,
    'Seq5': 1,
    'Seq6' :1
}

# Convert the alignment dictionary to a MultipleSeqAlignment object
alignment = MultipleSeqAlignment([
    SeqRecord(Seq(seq), id=name)
    for name, seq in alignment_data.items()
])

# Calculate the distance matrix
calculator = DistanceCalculator()
distance_matrix = calculator.get_distance(alignment)

# Construct the UPGMA tree
constructor = DistanceTreeConstructor()
tree = constructor.upgma(distance_matrix)

# Assign colors to clusters
cluster_colors = {
    1: 'red',
    2: 'blue',
    3: 'green'
}

# Customize the tree visualization
for terminal in tree.get_terminals():
    cluster = cluster_assignments.get(terminal.name)
    if cluster:
        terminal.color = cluster_colors.get(cluster)

# Print and visualize the tree
Phylo.draw(tree)

#To handle sequences of different lengths


from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Example sequence alignment
alignment_data = {
    'Seq1': 'AGCCG',
    'Seq2': 'AGAACG',
    'Seq3': 'ATAAAAAAG',
    'Seq4': 'ATAAAAAACCCGAG',
    'Seq5': 'ACGGGAG',
    'Seq6': 'ACTGGGGGTTTTTGAG'
}

# Cluster assignments for each sequence
cluster_assignments = {
    'Seq1': 3,
    'Seq2': 3,
    'Seq3': 2,
    'Seq4': 2,
    'Seq5': 4,
    'Seq6': 1
}

# Find the maximum length of sequences
max_length = max(len(seq) for seq in alignment_data.values())

# Pad sequences to the maximum length
alignment_data_padded = {
    name: seq.ljust(max_length, '-')
    for name, seq in alignment_data.items()
}

# Convert the padded alignment dictionary to a MultipleSeqAlignment object
alignment = MultipleSeqAlignment([
    SeqRecord(Seq(seq), id=name)
    for name, seq in alignment_data_padded.items()
])

# Calculate the distance matrix
calculator = DistanceCalculator()
distance_matrix = calculator.get_distance(alignment)

# Construct the UPGMA tree
constructor = DistanceTreeConstructor()
tree = constructor.upgma(distance_matrix)

# Assign colors to clusters
cluster_colors = {
    1: 'red',
    2: 'blue',
     3: 'green',
    4: 'yellow'
}

# Customize the tree visualization
for terminal in tree.get_terminals():
    cluster = cluster_assignments.get(terminal.name)
    if cluster:
        terminal.color = cluster_colors.get(cluster)

# Print and visualize the tree
Phylo.draw(tree)
