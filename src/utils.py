from Bio.SeqFeature import FeatureLocation
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import os

Entrez.email = "mokrane.mokhtari.etu@univ-lille.fr"

# Fonction pour extraire les positions des CDS à partir d'un enregistrement GenBank
def find_cds(seq_record):
    cds_positions = []
    for feature in seq_record.features:
        if feature.type == "CDS":  # Vérifier si c'est un CDS
            start = int(feature.location.start)
            end = int(feature.location.end)
            cds_positions.append((start, end))
    return cds_positions

# Fonction pour récupérer l'identifiant du gène à partir de l'ARNm en utilisant elink
def mrna_to_gene(acc):
    """
    Convertit un identifiant d'ARNm en identifiant de gène via NCBI Entrez.
    """
    try:
        handle = Entrez.esearch(db="gene", term=acc)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            raise ValueError(f"Aucun gène correspondant à l'ARNm {acc} n'a été trouvé.")

        gene_id = record["IdList"][0]
        print(f"🔍 ID du gène récupéré pour {acc} : {gene_id}")  # Debug
        return gene_id

    except Exception as e:
        print(f"❌ Erreur lors de la récupération du gène pour l'ARNm {acc}: {str(e)}")
        return None


# Fonction pour récupérer des informations génomiques sur un gène via esummary
def get_genomic_info_from_gene(gene_id):
    try:
        handle = Entrez.esummary(db="gene", id=mrna_to_gene(gene_id))
        record = Entrez.read(handle)
        handle.close()

        genomic_info = record["DocumentSummarySet"]["DocumentSummary"][0]["GenomicInfo"][0]
        chrom_accession = genomic_info["ChrAccVer"]
        start_pos = int(genomic_info["ChrStart"])
        end_pos = int(genomic_info["ChrStop"])
        strand = 1 if start_pos < end_pos else -1  # Déduction du brin

        return chrom_accession, start_pos, end_pos, strand
    except Exception as e:
        print(f"Erreur lors de la récupération des informations génomiques pour le gène {gene_id} : {str(e)}")
        return None, None, None, None

# Fonction pour récupérer la séquence amont du gène

def upstream_gene_seq(gene_id, length):
  
    
    handle = Entrez.esummary(db="gene", id=mrna_to_gene(gene_id))
    result = Entrez.read(handle)
    handle.close()
    
    genomic_info = result["DocumentSummarySet"]["DocumentSummary"][0].get("GenomicInfo", [])
    
    if not genomic_info:
        print(f"Aucune information génomique trouvée pour le gène {gene_id}.")
        return None

    chrom_accession = genomic_info[0]["ChrAccVer"]  # ID du chromosome
    start = int(genomic_info[0]["ChrStart"])        # Position de début du gène
    end = int(genomic_info[0]["ChrStop"])           # Position de fin du gène
    strand = int(genomic_info[0].get("Strand", 1))  # Sens du gène (1 = direct, -1 = inverse)

    print(f"Accession : {chrom_accession}, Start : {start}, End : {end}, Strand : {strand}")

    if start < end :
        endPos = start - length +1
        handle1 = Entrez.efetch(db="nucleotide", id = chrom_accession, rettype = "gb", retmode ="texte", seq_start=str(start), seq_stop=str(endPos), strand ="1")

    else: 
        endPos = start + length +1
        handle1 = Entrez.efetch(db="nucleotide", id = chrom_accession, rettype = "gb", retmode ="texte", seq_start=str(start+2), seq_stop=str(endPos), strand ="2")
    

    record1 = SeqIO.read(handle1, 'gb')
    sequence = record1.seq





    return sequence

def frequency_to_pssm(frequency_matrix):
  
    pwm = frequency_matrix.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()
    
    return pssm


def search_pssm_in_sequence(pssm, sequence):
    seq_len = len(sequence)
    motif_len = len(pssm)
    occurrences = []
    
    # Parcours de la séquence avec la PSSM
    for i in range(seq_len - motif_len + 1):
        subseq = sequence[i:i+motif_len]
        score = 0
        
        # Calcul du score de la sous-séquence par rapport à la PSSM
        for j in range(motif_len):
            base = subseq[j]
            if base == 'A':
                score += pssm[0, j]
            elif base == 'C':
                score += pssm[1, j]
            elif base == 'G':
                score += pssm[2, j]
            elif base == 'T':
                score += pssm[3, j]
        
        # Définir un seuil pour l'occurrence
        if score > 0:  # Ajustez ce seuil selon vos besoins
            occurrences.append(i)
    
    return occurrences

def download_promotors(mrna_ids, promoter_size, output_dir="data"):

    os.makedirs(output_dir, exist_ok=True)

    for mrna_id in mrna_ids:
        print(f"Téléchargement du promoteur pour {mrna_id}...")
        try:
            handle = Entrez.efetch(db="nucleotide", id=mrna_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            promoter_seq = record.seq[:promoter_size]

            filename = os.path.join(output_dir, f"{mrna_id}_{promoter_size}.fa")
            with open(filename, "w") as fasta_file:
                fasta_file.write(f">{mrna_id}_promoter_{promoter_size}\n{promoter_seq}\n")

            print(f"Enregistré dans {filename}")
        except Exception as e:
            print(f"Erreur avec {mrna_id}: {e}")

