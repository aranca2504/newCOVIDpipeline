import Bio
from Bio import SeqIO
import fastaparser
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

for record in SeqIO.parse("sequence.fasta", "fasta"):
    print(record.id)

### Muestra de referencia
with open("sequence.fasta") as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            # seq is a FastaSequence object
            print('ID:', seq.id)
            print('Descripción:', seq.description)
            print('Secuencia:', seq.sequence_as_string())
            #Filtrar las secuencias para asegurarnos de que realmente son completas
            print()

#Alinearlas utilizando MAFFT
mafft_cline = MafftCommandline(input="sequence.fasta")
stdout, stderr = mafft_cline()
with open("aligned.fasta","w") as handle:
     handle.write(stdout)
#Almacenar el resultado
alignment=AlignIO.read("sequence.fasta","fasta")
print(alignment)