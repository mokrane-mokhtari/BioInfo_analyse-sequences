import argparse

def main():
    parser = argparse.ArgumentParser(description="Script pour tester les arguments de ligne de commande.")

    parser.add_argument('-t', '--threshold', type=float, help="Seuil de score pour la matrice", required=True)
    
    parser.add_argument('-l', '--promotor-length', type=int, default=1000, help="Longueur du promoteur (par défaut 1000)")
    
    parser.add_argument('-w', '--window-size', type=int, default=40, help="Longueur de la fenêtre glissante (par défaut 40)")
    
    parser.add_argument('mRNA_ids', nargs='+', help="Liste des identifiants Genbank des mRNA")

    args = parser.parse_args()

    print(f"Seuil de score: {args.threshold}")
    print(f"Longueur du promoteur: {args.promotor_length}")
    print(f"Taille de la fenêtre glissante: {args.window_size}")
    print(f"Identifiants mRNA: {args.mRNA_ids}")

if __name__ == "__main__":
    main()
