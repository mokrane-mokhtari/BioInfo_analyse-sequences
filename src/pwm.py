from Bio import motifs
from Bio.Seq import Seq
from Bio import Entrez
from utils import mrna_to_gene, upstream_gene_seq


def pwm2pssm(freq_matrix, pseudocount=0.1):
    pwm = freq_matrix.normalize(pseudocount)
    pssm = pwm.log_odds()
    return pssm

def scan_sequence(pssm, sequence, threshold):
    occurrences = []
    length = len(sequence)
    for position, score in pssm.search(sequence, threshold=threshold, both=True):
        aposition = position - length
        if aposition >= -length:
            occurrences.append((aposition, score))
    return occurrences


def scan_all_sequences(pssm, sequences, threshold):
   
    scan_results = []
    for seq in sequences:
        scan_results.append(scan_sequence(pssm, seq, threshold))
    return scan_results


def score_window(scan_results, start, end):
   
    total_score = 0
    
    for occurrences in scan_results:
        for pos, score in occurrences:
            if start <= pos <= end:
                poids = abs(score)  
                total_score += score * poids
    
    return total_score



def best_window(scan_results, window_size):
    best_start, best_end = None, None
    best_score = float('-inf')

    all_occurrences = []
    
    for occurrences in scan_results:
        all_occurrences.extend([pos for pos, _ in occurrences])

    if not all_occurrences:
        return None, None, None
    all_occurrences = sorted(set(all_occurrences))


  
    for occ in all_occurrences:
        start = occ - window_size // 2  
        end = occ + window_size // 2
        current_score = score_window(scan_results, start, end)      

        print(f"fenêtre autour de l'occurrence à  {occ}: Début = {start}, Fin = {end}, Score = {current_score}")
        if current_score > best_score:
            best_score = current_score
            best_start, best_end = start, end


    
    return best_start, best_end, best_score

if __name__ == "__main__":
    Entrez.email = "mokrane.mokhtari.etu@univ-lille.fr"
    length = 500
    liste = ["NM_000451"]
    l = []
    file = "../data/MA0083.3.jaspar"
    threshold = -20
   

    for id in liste: 
        gene_id = mrna_to_gene(id)
        seq = upstream_gene_seq(gene_id, length)
        seq.id = id
        l.append(seq)

    with open(file) as handle:
        matrix = motifs.read(handle, "jaspar")
        pssm = pwm2pssm(matrix.counts, 0.1)

        results = scan_all_sequences(pssm, l, threshold)
        print("scan des sequences :", results)

        # Calculer le score pour une fenêtre spécifiée
        score_result = score_window(results, -200, -100)  
        print(f"Score de la fenêtre: {score_result}")
        
        # Recherche de la meilleure fenêtre
        best_window_result = best_window(results, window_size=40)  
        print("La meilleure fenêtre trouvée  :", best_window_result)

    
    