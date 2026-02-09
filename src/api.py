import os
from Bio import Entrez, SeqIO, motifs
from utils import frequency_to_pssm
from utils import find_cds, mrna_to_gene, upstream_gene_seq
from utils import search_pssm_in_sequence
from pwm import pwm2pssm , scan_sequence

# email 
Entrez.email = "mokrane.mokhtari.etu@univ-lille.fr"

# Numéro d'accession à récupérer
accession = "NM_007389"

# Partie 1 : Récupération de la séquence avec efetch
def fetch_sequence(accession, rettype="gb"):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype=rettype, retmode="text")
    record = SeqIO.read(handle, rettype)
    return record

# Récupération des fichiers
genbank_record = fetch_sequence(accession, "gb")
fasta_record = fetch_sequence(accession, "fasta")

# Affichage des résultats
print("ID de la séquence GenBank :", genbank_record.id)
print("Longueur de la séquence :", len(genbank_record.seq))

print("\nID de la séquence FASTA :", fasta_record.id)
print("Longueur de la séquence :", len(fasta_record.seq))


output_dir = "../data2"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Le dossier '{output_dir}' a ete cree.")

# Sauvegarde les fichiers localement
SeqIO.write(genbank_record, os.path.join(output_dir, "NM_007389.gb"), "genbank")
SeqIO.write(fasta_record, os.path.join(output_dir, "NM_007389.fasta"), "fasta")

print("\nLes fichiers ont été sauvegardés dans le dossier 'data2/'.")

cds_from_fetch = find_cds(genbank_record)
print("\nCDS extraite du fichier GenBank récupéré :", cds_from_fetch)

data_dir = "../data"
genbank_file_manual = os.path.join(data_dir, "NM_007389.gb")
local_genbank_record = SeqIO.read(genbank_file_manual, "genbank")

cds_from_local = find_cds(local_genbank_record)

if cds_from_fetch == cds_from_local:
    print("\nLa CDS extraite du fichier GenBank récupéré est identique à celle du fichier GenBank local.")
else:
    print("\nLa CDS extraite du fichier GenBank récupéré n'est PAS identique à celle du fichier GenBank local.")

# Récupérer l'ID du gène
try:
    gene_id = mrna_to_gene(accession)
    print(f"L'identifiant du gène correspondant à l'ARNm {accession} est {gene_id}.")
except ValueError as e:
    print(f"Erreur : {e}")


length = 1000  
upstream_sequence = upstream_gene_seq(gene_id, length)

if upstream_sequence:
    print(f"Séquence amont du gène {gene_id} :\n{upstream_sequence}")
else:
    print("Impossible de récupérer la séquence amont.")


## dernier partie 
fichier_jaspar = "../data/JASPAR2024.txt"

# Lecture de toutes les matrices depuis le fichiers JASPAR
with open(fichier_jaspar) as handle:
    motifs_list = list(motifs.parse(handle, "jaspar"))

# Verification du nombre de matrices lues
print(f" Nombre de matrices lues : {len(motifs_list)}")


if motifs_list:
    motif = motifs_list[0]  
    
    # Obtenir la matrice de type FrequencyPositionMatrix
    print(f" Motif ID : {motif.name}")
    print(f" Matrice de fréquences (FPM) :\n{motif.counts}")
    
    # Transformer en PWM 
    pwm = motif.counts.normalize(pseudocounts=0.5)
    print(f" Position-Weight Matrix (PWM) :\n{pwm}")

    # Transformer en PSSM 
    pssm = pwm.log_odds()
    print(f" Position-Specific Scoring Matrix (PSSM) :\n{pssm}")
    frequency_matrix = motif.counts
    
    # Convertir la matrice de fréquence en PSSM
    pssm = frequency_to_pssm(frequency_matrix)
    
    # Afficher la matrice PSSM
    print(f" Position-Specific Scoring Matrix (PSSM) en passant drct par la foncton:\n{pssm}")

    #appel pour rechercher les occurences d'une PSSM dans une séquence 
    sequence_to_search = fasta_record.seq  # on peut mettre nimporte quelle sequence quon veut 
    occurrences = search_pssm_in_sequence(pssm, sequence_to_search)

    # Affichage des occurrences
    print(f"Les occurrences du motif dans la séquence sont aux positions : {occurrences}")

    #on  transfome en PSSM en utilsant pwm2pssm 
    pssm_pwm2pssm = pwm2pssm(frequency_matrix)
    print(f" Position-Specific Scoring Matrix (PSSM) en utilisant pwm2pssm :\n{pssm_pwm2pssm}")

    #pour appeler la fonction scan_sequence
    seuil = -2.0
    results = scan_sequence(pssm_pwm2pssm,sequence_to_search,seuil)
    print(f"Les occurrences du motif dans la séquence sont aux positions et scores : {results}")
    