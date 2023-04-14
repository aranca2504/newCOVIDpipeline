import os
import subprocess
from io import StringIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline
import shutil
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
    # Inicializamos listas donde almacenaremos los numeros de cambios por posicion
    nogap_num_changes = []
    absolute_num_changes = []
    # Iterar por cada columna del alineamiento
    for i in tqdm(range(alignment.get_alignment_length())):
        column = alignment[:, i]  # nn-ntnn

        # Generar un diccionario con las entradas y sus ocurrencias
        column_counts = {}  # {'n': 3, '-': 2, 't': 1}
        for char in column:
            if char in column_counts:
                column_counts[char] += 1
            else:
                column_counts[char] = 1

        # CONTAR DIFERENCIAS EN CADA COLUMNA
        # con gaps
        max_val = max(column_counts.values())  # obtenemos el valor máximo del diccionario
        absolute_num_changes_calc = sum([v for k, v in column_counts.items() if v != max_val])
        absolute_num_changes.append(absolute_num_changes_calc)
        # sin gaps
        nogap_column_counts = {key: value for key, value in column_counts.items() if key != '-'}
        nogap_num_changes_calc = sum([v for k, v in nogap_column_counts.items() if v != max_val])
        nogap_num_changes.append(nogap_num_changes_calc)
        # Calcular el número de posiciones que tienen algún gap en el alineamiento

    # PLOT NUMERO DE CAMBIOS POR POSICION
    # ignorar el tipico warning de matplotlib:
    warnings.filterwarnings("ignore", message="Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.")
    posiciones = range(len(nogap_num_changes))

    # Lo mismo pero teniendo en cuenta los gaps
    warnings.filterwarnings("ignore", message="Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.")
    posiciones_abs = range(len(absolute_num_changes))

    print("GRAFICANDO DISTRIBUCIONES...")
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(8, 8))
    # Trazar la primera gráfica en el primer subplot
    ax1.bar(posiciones_abs, absolute_num_changes)
    ax1.set_title('Número de cambios con gap:' + str(sum(absolute_num_changes)))
    ax1.set_xlabel('Posición')
    ax1.set_ylabel('nº de cambios')
    # Trazar la segunda gráfica en el segundo subplot
    ax2.bar(posiciones, nogap_num_changes)
    ax2.set_title('Número de cambios sin gap:' + str(sum(nogap_num_changes)))
    ax2.set_xlabel('Posición')
    ax2.set_ylabel('nº de cambios')
    # Mostrar gráficas
    plt.show()

def checkAlign(alignment, new_aligned_fasta):
    """
    A function to check the quality of a multiple sequence alignment and remove highly gapped sequences.
    :param alignment: A Biopython alignment object containing the multiple sequence alignment.
    :param new_aligned_fasta: The path to the output FASTA file for the new alignment.
    """

    plotChanges(alignment)
    seq_index_to_remove = []
    print("Comprobando delecciones")
    for j in range(alignment.__len__()):

        fila = alignment[j, :]
        num_guiones = fila.count('-')
        length = fila.__len__()
        prop_delec = num_guiones / length

        # ELIMINAR SECUENCIAS CON MUCHOS GAPS
        if prop_delec > 0.05:
            seq_index_to_remove.append(j)  # Por ejemplo, eliminar la última secuencia
            print("Secuencia eliminada:" + str(j))
        new_align = MultipleSeqAlignment([])
    for i, record in enumerate(alignment):
        if not i in seq_index_to_remove:
            new_align.append(record)
    # Guardar la nueva alineación en un archivo FASTA
    AlignIO.write(new_align, new_aligned_fasta, "fasta")

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
        # Verifica que el nombre del archivo comience con "new_aligned_sequences.fasta."
        if file_name.startswith(new_aligned_fasta + "."):
            # Construye la ruta completa del archivo original
            file_path = os.path.join(source_dir, file_name)
            # Construye la ruta completa donde se moverá el archivo
            dest_path = os.path.join(output_iqtree_dir, file_name)
            # Mueve el archivo al nuevo directorio
            shutil.move(file_path, dest_path)

    return tree

if __name__ == '__main__':
    #DEFINE DIRECTORYS
    input_dir = "input_sequences"  # Ruta al directorio con todas las secuencias fasta en diferentes archivos
    output_sequences = "output_sequences" # Ruta al directorio con las secuencias de salida (conc_fasta, aligned_fasta, new_aligned_fasta)
    output_iqtree_dir = "output_iqtree"
    path_iqTree = "/Users/macbook/Downloads/iqtree-2.2.0-MacOSX/bin/iqtree2"  # Ruta al ejecutable de igTree


    #DEFINE FILE NAMES
    conc_fasta = "conc_sequences.fasta"  # Ruta del archivo que se creará con las secuencias de input concatenadas
    aligned_fasta = "aligned_sequences.fasta"  # Ruta del archivo con las secuencias alineadas con MAFFT
    new_aligned_fasta = "new_aligned_sequences.fasta"  # Ruta del archivo con las secuencias realineadas
    os.makedirs(output_sequences, exist_ok=True)
    conc_fasta = os.path.join(output_sequences,conc_fasta)
    aligned_fasta = os.path.join(output_sequences,aligned_fasta)
    new_aligned_fasta = os.path.join(output_sequences,new_aligned_fasta)

    #RUN FUNCTIONS
    concatenate_fasta_files(input_dir, conc_fasta)
    alignment, stderr = mafft_align_multiple(conc_fasta, aligned_fasta)
    checkAlign(alignment, new_aligned_fasta)
    tree = createTree(path_iqTree, new_aligned_fasta,output_sequences, output_iqtree_dir)

    #MANAGE TREE
    tree.render("arbol.png", w=400)
    print(stderr)
    print(tree)
    tree.show()
