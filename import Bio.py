import Bio
from Bio import Entrez
from Bio import SeqIO
import fastaparser

for record in SeqIO.parse("sequence.fasta", "fasta"):
    print(record.id)

### Muestra de Wuhan
with open("sequence.fasta") as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            # seq is a FastaSequence object
            print('ID:', seq.id)
            print('Descripción:', seq.description)
            print('Secuencia:', seq.sequence_as_string())
            print()