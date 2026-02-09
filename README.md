# Projet Bioinformatique
## Realisé par :
- MOKHTARI Mokrane
- BENANI Dalya

Ce projet est une série de scripts Python destinés à l'analyse de séquences génétiques, en particulier pour l'identification de motifs de liaison de facteurs de transcription (TFBS) dans des séquences d'ADN. Le projet inclut plusieurs modules, chacun ayant des fonctionnalités spécifiques pour le traitement et l'analyse des données biologiques.

## Description générale

Le projet est principalement composé de scripts qui exécutent des analyses sur des séquences génétiques. Les principaux objectifs sont :

- Extraction de séquences génétiques à partir de fichiers au format GenBank et FASTA.
- Analyse de motifs de liaison de facteurs de transcription (TFBS) en utilisant des matrices de position de motifs (PSSM) provenant de bases de données comme JASPAR.
- Tests et affichage des résultats sous forme de fichiers JSON pour une analyse détaillée des occurrences et des scores.

## Structure du projet

- Dossier data/ : Ce dossier contient les données nécessaires pour les analyses, comme des fichiers au format .fa (FASTA), .txt, .gb (GenBank), etc. Ces fichiers servent de base pour le traitement des séquences génétiques.

- Dossier src/ : Ce dossier contient les scripts Python qui réalisent les différentes étapes d'analyse des séquences génétiques :

- api.py : Un script pour récupérer et afficher des données génétiques à partir de fichiers comme GenBank ou FASTA.

- scan_pwm.py : Un script qui scanne des motifs de liaison de facteurs de transcription (TFBS) dans des séquences génétiques.

- pwm.py : Un script pour effectuer des calculs sur les matrices de position des motifs (PWM).

- test_args.py : Un script pour tester les arguments passés en ligne de commande, comme les seuils et autres paramètres pour les analyses.

- putativeTFBS.py : Un script qui prédit les TFBS dans les séquences génétiques en scannant avec des matrices de motifs.

- resultats.json : Ce fichier contient les résultats des analyses sous format JSON. Il inclut les fenêtres où des facteurs de transcription ont été détectés, avec des informations sur leur score et leur position.

- README.md : Ce fichier qui contient la documentation du projet, expliquant son fonctionnement, son installation, et son utilisation.

## Utilisation
# Exécution de l'API 

Pour récupérer les données de séquences génétiques depuis un fichier GenBank ou FASTA, exécutez le script api.py : 

```
python3 api.py

``` 

# Exécution du script scan_pwm.py : 

Ce script permet de scanner un fichier de motifs de liaison de facteurs de transcription (PSSM) pour une séquence donnée. Voici un exemple de commande pour l'exécuter : 

```
python3 scan_pwm.py ../data/JASPAR2024.txt ../data/NM_007389 1000 -2.0
```
Cette commande prend un fichier de motifs PSSM (JASPAR2024.txt) et une séquence d'ADN (NM_007389), puis recherche les motifs dans la séquence en utilisant les paramètres donnés.

# Exécution du script pwm.py: 

Le script pwm.py effectue des calculs liés aux matrices de motifs de liaison. Vous pouvez l'exécuter de cette manière :

```
python3 pwm.py ../data/NM_007389.fasta
```

# Exécution du script test_args.py: 

Ce script est utilisé pour tester les arguments en ligne de commande, en particulier les valeurs pour le seuil de score, la longueur de la séquence et la taille de la fenêtre. Voici la commande pour l'exécuter :

```
python3 test_args.py -t 5.5 -l 1200 -w 50 NM_001 NM_002 NM_003
```
# Exécution du script putativeTFBS.py : 

Ce script permet de prédire des TFBS en scannant des séquences génétiques avec une matrice de motifs. Pour exécuter ce script, utilisez la commande suivante : 

```
python3 putativeTFBS.py -m ../data/MA0083.3.jaspar -t -20 -s -10 NM_000451

```

Cette commande prend un fichier de matrice de motifs JASPAR (MA0083.3.jaspar), une séquence d'ADN (NM_000451 ), et divers paramètres pour prédire les TFBS, puis enregistre les résultats dans un fichier JSON (resultats.json).

## Format des résultats (JSON): 
Ce format JSON représente les résultats d'une analyse de motifs de liaison de facteurs de transcription (TFBS) dans une séquence génétique : 


- window_num : Le numéro de la fenêtre analysée .

- tf_name : Le nom du facteur de transcription détecté.

- window_range : La plage de positions de la fenêtre dans la séquence.

- window_score : Le score de la fenêtre indiquant la force de la correspondance du motif .

- num_occurrences : Le nombre d'occurrences du motif détectées dans la fenêtre .

- occurrences : Détails de chaque occurrence trouvée, incluant :

- mrna_id : Identifiant du mRNA .

- position : Position de l'occurrence dans la séquence .

- score : Le score de l'occurrence .
