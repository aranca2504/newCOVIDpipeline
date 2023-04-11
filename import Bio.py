import Bio
from Bio import SeqIO
import fastaparser
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import os
from io import StringIO
import matplotlib.pyplot as plt

for record in SeqIO.parse("sequence.fasta", "fasta"):
    print(record.id)

### PASO 1. Muestra de referencia
with open("sequence.fasta") as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            # seq is a FastaSequence object
            print('ID:', seq.id)
            print('Descripción:', seq.description)
            print('Secuencia:', seq.sequence_as_string())
            #Filtrar las secuencias para asegurarnos de que realmente son completas
            print()

#PASO 2. Alinearlas utilizando MAFFT
mafft_cline = MafftCommandline(input="sequence.fasta")
stdout, stderr = mafft_cline()
with open("aligned.fasta","w") as handle:
     handle.write(stdout)
#Almacenar el resultado
alignment=AlignIO.read("sequence.fasta","fasta")
print(alignment)

#PASO 3.Comprobar alineamiento
# Inicializar el diccionario de cambios por columna
changes_per_column = {}

# Iterar por cada columna del alineamiento
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]

    # Contar los diferentes caracteres (incluyendo gaps) en la columna
    column_counts = {}
    for char in column:
        if char in column_counts:
            column_counts[char] += 1
        else:
            column_counts[char] = 1

    # Contar el número total de cambios (excluyendo gaps) en la columna
    num_changes = sum(count for char, count in column_counts.items() if char != "-") - 1

    # Agregar el número de cambios a la lista correspondiente en el diccionario
    if num_changes in changes_per_column:
        changes_per_column[num_changes] += 1
    else:
        changes_per_column[num_changes] = 1

        x_vals = list(changes_per_column.keys())
        y_vals = list(changes_per_column.values())

        plt.bar(x_vals, y_vals)
        plt.xlabel("Número de cambios por columna")
        plt.ylabel("Número de posiciones")
        plt.show()

        # Calcular el número de posiciones que tienen algún gap en el alineamiento
        num_gap_positions = 0
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            if "-" in column:
                num_gap_positions += 1

                print("El número de posiciones con gaps es:", num_gap_positions)

                # Si hay secuencias concretas que tienen muchos gaps, eliminarlas
                seq_index_to_remove = alignment.get_alignment_length()  # Por ejemplo, eliminar la última secuencia
                new_alignment = alignment[:, :seq_index_to_remove] + alignment[:, seq_index_to_remove+1:]

                # Guardar la nueva alineación en un archivo FASTA
                AlignIO.write(new_alignment, "new_aligned.fasta", "fasta")
