import os
import subprocess
from io import StringIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
import numpy as np
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
from tqdm import tqdm
import PyQt5 #necesario para ete3
import matplotlib.pyplot as plt
import warnings

def concatenate_fasta_files(input_dir, conc_fasta):
    """
    A function to concatenate multiple fasta files in a directory into a single output file.
    :param input_dir: The path to the directory containing the input fasta files.
    :param conc_fasta: The path to the output concatenated fasta file.
    """
    input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.fasta')]
    with open(conc_fasta, 'w') as out_file:
        for file in input_files:
            with open(file, 'r') as in_file:
                for line in in_file:
                    out_file.write(line)

def mafft_align_multiple(conc_fasta, aligned_fasta):
    """
    A function to align multiple fasta files using MAFFT.
    Prints the resultant alignment
    :param conc_fasta: The path to the concatenated fasta file
    :param aligned_fasta: The path to the output aligned fasta file.
    :return alignment: Resultant alignment class-MultipleSeqAlignment
    :return stderr: string - Historial of the strategy followed for the alignment
    """
    mafft_cline = MafftCommandline(input=conc_fasta)
    aligned_seqs, stderr = mafft_cline()  # stderr es ek

    # Parse the aligned sequences with AlignIO
    alignment = AlignIO.read(StringIO(aligned_seqs), "fasta")

    # Write the aligned sequences to a new file
    with open(aligned_fasta, "w") as handle:
        AlignIO.write(alignment, handle, "fasta")

    return alignment, stderr

def plotChanges(alignment):
    """
    A function to create a .png file in which are plotted the distribution of Divergence between sequences
    :param alignment: Resultant alignment class-MultipleSeqAlignment
    """

    # Initialize lists where we will store the number of changes per position
    absolute_char_changes = [] # GENERAL CHARACTER CHANGES (e.g. 'nn-nttn' = 2): max achievable value = 4
    nogap_char_changes = [] # CHARACTER CHANGES DISREGARDING GAPS (e.g. 'nn-nttn' = 1): max achievable value = 3
    absolute_num_changes = [] # NUMBER OF DISTINCT CHARS IN COLUMN FROM MAX (e.g. 'nn-nttn' = 3)
    nogap_num_changes = [] # NUMBER OF DISTINCT CHARS IN COLUMN FROM MAX DISREGARDING GAPS (e.g. 'nn-nttn' = 2)

    # Iterate over each column in the alignment
    for i in tqdm(range(alignment.get_alignment_length())):
        column = alignment[:, i]  # nn-ntnn

        # Generate a dictionary with the entries and their occurrences
        column_counts = {} # {'n': 3, '-': 2, 't': 1}
        for char in column:
            if char in column_counts:
                column_counts[char] += 1
            else:
                column_counts[char] = 1

        # COUNT DIFFERENCES IN EACH COLUMN
        # with gaps
        max_val = max(column_counts.values()) # get the maximum value from the dictionary
        absolute_char_changes_calc = sum([v for k, v in column_counts.items() if v != max_val]) # count differences
        absolute_char_changes.append(absolute_char_changes_calc) # store in the vector
        absolute_num_changes_calc = sum(column_counts.values()) - max_val # count differences (sum of all - most common)
        absolute_num_changes.append(absolute_num_changes_calc)
        # without gaps (same as with gaps except for the first line)
        nogap_column_counts = {key: value for key, value in column_counts.items() if key != '-'} # get column without gaps
        nogap_char_changes_calc = sum([v for k, v in nogap_column_counts.items() if v != max_val])
        nogap_char_changes.append(nogap_char_changes_calc)
        nogap_num_changes_calc = sum(nogap_column_counts.values()) - max_val
        nogap_num_changes.append(nogap_num_changes_calc)

    absolute_char_changes = np.array(absolute_char_changes) # convert to np.array
    nogap_char_changes = np.array(nogap_char_changes)

    # PLOT NUMBER OF CHANGES PER POSITION
    # ignore typical matplotlib warning:
    warnings.filterwarnings("ignore", message="Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.")
    positions = np.array(range(len(nogap_char_changes)))
    # same but with gaps
    warnings.filterwarnings("ignore", message="Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.")
    positions_abs = np.array(range(len(absolute_char_changes)))

    print("PLOTTING DISTRIBUTIONS...")
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 10))
    # plot the first graph in the first subplot
    axs[0, 0].bar(positions_abs, absolute_char_changes)
    axs[0, 0].set_title('Number of different BP with gaps:' + str(sum(absolute_char_changes)))
    axs[0, 0].set_xlabel('Position')
    axs[0, 0].set_ylabel('Number of changes')
    # Trazar la segunda gráfica en el segundo subplot
    axs[1, 0].bar(positions, nogap_char_changes)
    axs[1, 0].set_title('Number of different BP without gaps:' + str(sum(nogap_char_changes)))
    axs[1, 0].set_xlabel('Position')
    axs[1, 0].set_ylabel('Number of changes')
    # Trazar la tercera
    axs[0, 1].bar(positions_abs, absolute_num_changes)
    axs[0, 1].set_title('Number of total variations with gaps:' + str(sum(absolute_num_changes)))
    axs[0, 1].set_xlabel('Position')
    axs[0, 1].set_ylabel('Number of changes')
    # Trazar la cuarta
    axs[1, 1].bar(positions, nogap_num_changes)
    axs[1, 1].set_title('Number of total variations without gap:' + str(sum(nogap_num_changes)))
    axs[1, 1].set_xlabel('Position')
    axs[1, 1].set_ylabel('Number of changes')
    # Mostrar gráficas
    plt.subplots_adjust(hspace=0.4, wspace=0.4)#ajustar espacio
    plt.savefig("distributions.png")

def checkAlign(alignment, new_aligned_fasta):
    """
    A function to check the quality of a multiple sequence alignment and remove highly gapped sequences.
        1) Erases sequences with excesive gaps
        2) Removes columns with gaps
        3) Creates a new file containing the new alignment.
    :param alignment: A Biopython alignment object containing the multiple sequence alignment.
    :param new_aligned_fasta: The path to the output FASTA file for the new alignment.
    :return new_align_1: (MultipleSeqAlignment) The alignment with correction 1)
    :return new_align_2: (MultipleSeqAlignment) The alignment with corrections 1) and 2)
    """

    # 1)Erases sequences with excesive gaps
    ids_og = []
    seq_index_to_remove = []
    for record in alignment:
        ids_og.append(record.id)

    print("Checking sequences with a high number of detetions...")
    for j in range(alignment.__len__()):
        fila = alignment[j, :]
        num_guiones = fila.count('-')
        length = fila.__len__()
        prop_delec = num_guiones / length

        # REMOVE SEQUENCES WITH MANY GAPS
        if prop_delec > 0.05:
            seq_index_to_remove.append(j)
            print("Sequence removed:" + str(j))
        else:
            print("Sequence not removed -> ID: "+ids_og[j])
    new_align_1 = MultipleSeqAlignment([]) # Create de new alignment
    for i, record in enumerate(alignment):
        if not i in seq_index_to_remove:
            new_align_1.append(record)

    # 2) Removes columns with gaps
    # Initialize parameters
    ids = []  # Store the sequence ids, we will need them to reconstruct the alignment
    for record in new_align_1:
        ids.append(record.id)
    col_index_to_remove = []

    print("Checking deletions, step 2: remove columns with deletions")
    for j in range(new_align_1.get_alignment_length()):
        col = new_align_1[:, j]
        # STORE COLUMNS INDICES WITH GAPS
        if '-' in col:
            col_index_to_remove.append(j)

    # Transpose the matrix
    alignment_array = np.array(
        [list(rec.seq) for rec in new_align_1])  # Convert the new_align_1 matrix into a NumPy array
    alignment_transposed = alignment_array.T  # Transpose the alignment matrix
    alignment_transposed_records = [  # Convert the transposed NumPy array into a new MultipleSeqAlignment object
        SeqRecord(Seq("".join(row)))
        for i, row in enumerate(alignment_transposed)]
    new_align_1_T = MultipleSeqAlignment(alignment_transposed_records)
    # Generate a new alignment with only sequences without gaps
    new_align_2_T = MultipleSeqAlignment([])
    for i, record in tqdm(enumerate(new_align_1_T), desc="Total iterations: " + str(len(new_align_1_T)) + "  //  Current iteration"):
        if not i in col_index_to_remove:
            new_align_2_T.append(record)

    # Transpose the obtained alignment to return it to its original form
    alignment_array_2 = np.array([list(rec.seq) for rec in new_align_2_T])  # Convert the alignment into a NumPy array
    alignment_transposed_2 = alignment_array_2.T  # Transpose the matrix
    alignment_transposed_records_2 = [  # Convert the transposed NumPy matrix into a new MultipleSeqAlignment object
        SeqRecord(Seq("".join(row)), id=ids[i])
        for i, row in enumerate(alignment_transposed_2)]
    new_align_2 = MultipleSeqAlignment(alignment_transposed_records_2)

    #Guardar el alignment
    with open(new_aligned_fasta, "w") as handle:
        AlignIO.write(new_align_2, handle, "fasta")

    return new_align_1, new_align_2

def createTree(path_igTree, new_aligned_fasta,output_sequences, output_iqtree_dir):
    """
    A function to create a phylogenetic tree using IgTree software.
    :param path_igTree: The path to the IgTree executable file.
    :param new_aligned_fasta: The path to the input aligned fasta file.
    :param output_iqtree_dir: The path to the output directory.
    :return: tree: An ete3 tree object.
    """
    os.makedirs(output_iqtree_dir, exist_ok=True)
    subprocess.run([path_igTree, "-s", new_aligned_fasta, "-m", "GTR+I", "-bb", "1000"])
    tree_name = new_aligned_fasta + ".treefile"
    tree = Tree(tree_name)

    # MOVE NEW ARCHIVES TO THE OUTPUT DIRECTORY
    source_dir = output_sequences
    # remove the directory path from the file name
    new_aligned_fasta = os.path.basename(new_aligned_fasta)

    for file_name in os.listdir(source_dir):
        # Check that the file name starts with "new_aligned_sequences.fasta."
        if file_name.startswith(new_aligned_fasta + "."):
            # Build the complete path of the original file
            file_path = os.path.join(source_dir, file_name)
            # Build the complete path where the file will be moved
            dest_path = os.path.join(output_iqtree_dir, file_name)
            # Move the file to the new directory
            shutil.move(file_path, dest_path)

    return tree

if __name__ == '__main__':
    #DEFINE DIRECTORYS
    input_dir = "input_sequences"  # Path to the directory with all fasta sequences in separate files
    output_sequences = "output_sequences"  # Path to the directory with the output sequences (conc_fasta, aligned_fasta, new_aligned_fasta)
    output_iqtree_dir = "output_iqtree"
    path_iqTree = "/Users/macbook/Downloads/iqtree-2.2.0-MacOSX/bin/iqtree2"  # Path to the igTree executable

    #DEFINE FILE NAMES
    conc_fasta = "conc_sequences.fasta"  # Path to the file that will be created with concatenated input sequences
    aligned_fasta = "aligned_sequences.fasta"  # Path to the file with sequences aligned with MAFFT
    new_aligned_fasta = "new_aligned_sequences.fasta"  # Path to the file with the realigned sequences
    os.makedirs(output_sequences, exist_ok=True)
    conc_fasta = os.path.join(output_sequences, conc_fasta)
    aligned_fasta = os.path.join(output_sequences, aligned_fasta)
    new_aligned_fasta = os.path.join(output_sequences, new_aligned_fasta)

    #RUN FUNCTIONS
    concatenate_fasta_files(input_dir, conc_fasta)  # concatenate all fasta files in the input_dir
    alignment, stderr = mafft_align_multiple(conc_fasta, aligned_fasta)  # align the concatenated fasta file using MAFFT
    new_alignment_1, new_alignment_2 = checkAlign(alignment, new_aligned_fasta)  # realign sequences and create a new alignment file
    plotChanges(alignment)  # create a plot showing the alignment changes
    tree = createTree(path_iqTree, new_aligned_fasta, output_sequences, output_iqtree_dir)  # create a phylogenetic tree using IgTree


    #MANAGE TREE
    tree.render("arbol.png", w=400)  # render the tree to a png file
    print(stderr)  # print any errors from the MAFFT alignment
    print(tree)  # print the tree object
    tree.show()  # show the tree in the console or in a viewer
