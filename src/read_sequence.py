from Bio import SeqIO
from utils import find_cds
from utils import mrna_to_gene
from compareseq import compare_sequences

#pour les chemins des fichiers
fasta_file = "../data/NM_007389.fasta"
genbank_file = "../data/NM_007389.gb"


#pour charger les fichiers
fasta_record = SeqIO.read(fasta_file, "fasta")
genbank_record = SeqIO.read(genbank_file, "genbank")

#pour verifier les sequence
print("Les séquences sont identiques :", fasta_record.seq == genbank_record.seq)

#conversion en chaine de caracteres
fasta_str = str(fasta_record.seq)
genbank_str = str(genbank_record.seq)

print("Premiers nucléotides de la séquence FASTA :", fasta_str[:100])

reverse_complement = fasta_record.seq.reverse_complement()
print("Reverse complement :", reverse_complement[:100])

# Affichage des annotations (uniquement dans GenBank)
print("Annotations disponibles :", genbank_record.annotations.keys())

# Nombre de features et première feature
print("Nombre de features :", len(genbank_record.features))
if len(genbank_record.features) > 0:
    first_feature = genbank_record.features[0]
    print("Première feature - Début :", first_feature.location.start, "Fin :", first_feature.location.end)

print(compare_sequences(fasta_str,genbank_str))

# Trouver les CDS
cds_list = find_cds(genbank_record)
print("CDS trouvés :", cds_list)