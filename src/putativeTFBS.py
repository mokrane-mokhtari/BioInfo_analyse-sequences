import argparse
import json
import os
from Bio import SeqIO, motifs
from pwm import pwm2pssm, scan_all_sequences
from utils import *  
def main():
    parser = argparse.ArgumentParser(description="Recherche d’occurrences de matrices PFM dans un ensemble de promoteurs.")
    
    parser.add_argument('-m', '--pfm', required=True, help="Fichier contenant la matrice JASPAR")
    parser.add_argument('-t', '--threshold', required=True, type=float, help="Seuil de score de matrice")
    parser.add_argument('-l', '--promotor-length', default=500, type=int, help="Longueur du promoteur (par défaut 500)")
    parser.add_argument('-w', '--window-size', default=40, type=int, help="Longueur de la fenêtre glissante (par défaut 40)")
    parser.add_argument('-s', '--window-threshold', required=True, type=float, help="Seuil de score de fenêtre")
    parser.add_argument('-p', '--pseudocount', default=0.1, type=float, help="Valeur du pseudo-poids (par défaut 0.1)")
    parser.add_argument("mrna_ids", nargs="+", help="Liste d'identifiants mRNA.")

    args = parser.parse_args()

    with open(args.pfm) as handle:
        matrix = motifs.parse(handle, "jaspar")[0]

    pssm = pwm2pssm(matrix.counts, args.pseudocount)

    output_dir = "data"  
    download_promotors(args.mrna_ids, args.promotor_length, output_dir)

    sequences = []
    for mrna_id in args.mrna_ids:
        fasta_file = os.path.join(output_dir, f"{mrna_id}_{args.promotor_length}.fa")

        if os.path.exists(fasta_file):
            seq_record = SeqIO.read(fasta_file, "fasta")
            sequences.append(seq_record.seq)
        else:
            print(f"Impossible de trouver la séquence pour {mrna_id}.")

    scan_results = scan_all_sequences(pssm, sequences, threshold=args.threshold)
    window_results = []

    for i, seq in enumerate(sequences):
        mrna_id = args.mrna_ids[i]

        if not scan_results[i]:  
            print(f"Aucune occurrence trouvée pour {mrna_id}.")
            continue

        for pos, score in scan_results[i]:
            start = pos - (args.window_size // 2)
            end = pos + (args.window_size // 2)

            
            start = max(start, 0)  
            end = min(end, len(seq))  

            result = {
                "window_num": len(window_results) + 1,
                "tf_name": matrix.name,
                "window_range": f"[{start}:{end}]",
                "window_score": float(score),  
                "num_occurrences": 1,
                "occurrences": [{
                    "occurrence_num": 1,
                    "mrna_id": mrna_id,
                    "tf_name": matrix.name,
                    "position": int(pos),
                    "score": float(score)
                }]
            }
            
            window_results.append(result)

    for result in window_results:
        window_score = sum(occurrence['score'] for occurrence in result['occurrences']) / len(result['occurrences']) if result['occurrences'] else 0
        result['window_score'] = float(window_score)

    for result in window_results:
        print(f"{result['window_num']} {result['tf_name']} {result['window_range']} {result['window_score']:.2f} {result['num_occurrences']}")
        for occurrence in result["occurrences"]:
            print(f"{occurrence['occurrence_num']} {occurrence['mrna_id']} {occurrence['tf_name']} {occurrence['position']} {occurrence['score']:.2f}")
        print()

    output_file = "resultats.json"
    with open(output_file, "w") as json_file:
        json.dump(window_results, json_file, indent=4)

    print(f"Les résultats ont été enregistrés dans {output_file}.")

if __name__ == "__main__":
    main()
