import sys
from pwm import pwm2pssm, scan_sequence  # Import des fonctions nécessaires
from Bio import SeqIO
from Bio import motifs
from utils import upstream_gene_seq  # Suppression de search_pssm_in_sequence s'il n'est pas utilisé

# Vérifier qu'il y a bien 4 arguments passés en ligne de commande
if len(sys.argv) != 5:
    print("Erreur : Veuillez fournir 4 arguments :")
    print("Usage: python scan_pwm.py <fichier_matrice_jaspar> <identifiant_genbank> <taille_promotrice> <seuil_score>")
    sys.exit(1)

# Récupérer les arguments
fichier_matrice_jaspar = sys.argv[1]  # Fichier contenant la matrice JASPAR
identifiant_genbank = sys.argv[2]  # ID du gène ou séquence
taille_promotrice = int(sys.argv[3])  # Longueur de la région promotrice
seuil_score = float(sys.argv[4])  # Seuil pour considérer une occurrence

with open(fichier_matrice_jaspar, "r") as handle:
    motifs_list = list(motifs.parse(handle, "jaspar"))

if not motifs_list:
    print("Erreur : Aucun motif JASPAR chargé depuis le fichier fourni.")
    sys.exit(1)

for motif in motifs_list:
    print(f"Motif chargé : {motif.name}")

sequence = upstream_gene_seq(identifiant_genbank, taille_promotrice)
if sequence is None:
    print(f"Erreur : Impossible de récupérer la séquence promotrice pour {identifiant_genbank}.")
    sys.exit(1)

PSEUDO_POIDS = 0.1


for motif in motifs_list:
    pssm = pwm2pssm(motif.counts, PSEUDO_POIDS)  
    occurrences = scan_sequence(pssm, sequence, seuil_score)  

   
    for position, score in occurrences:     
        
        print(f"{motif.name} {position} {score}")

