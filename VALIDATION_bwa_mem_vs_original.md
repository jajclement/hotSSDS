# Validation : comparaison ancien pipeline (bwa aln/bwa-ra/sampe) vs nouveau (bwa mem)

## Objectif

Vérifier que la migration vers `bwa mem` + `bin/ITR_id_bwamem.pl` (voir [CHANGELOG_bwa_mem_migration.md](CHANGELOG_bwa_mem_migration.md)) reproduit fidèlement la classification ssDNA de l'ancien pipeline, en comparant les deux sur le même jeu de données.

## Méthode

Les deux pipelines ont été exécutés sur le même jeu de test (`tests/fastq/input.csv`, données HiSeq, puisque les données NovaSeq faisaient planter l'ancien pipeline, afin d'avoir une référence fonctionnelle des deux côtés) :

Mêmes paramètres (`conf/cluster.config`, `conf/mm10.json`, `--with_trimgalore`, profil `singularity`), répertoires de travail et de sortie séparés.

## Résultats

### 1. Volume de lectures en entrée de `parseITRs`

Les *unparsed sont les fichiers alignés bruts puis filtrés (secondaires, mauvaise qualité)

| Pipeline | Lignes (`*.unparsed.bed`) |
|---|---|
| Ancien | 13 820 |
| Nouveau | 21 989 (+59 %) |

Le nouveau pipeline conserve significativement plus de lectures jusqu'à l'étape de classification. 

Interprétation : `bwa mem` est plus sensible que `bwa aln`/`bwa-ra` (qui pouvait échouer silencieusement à aligner certaines paires avec ITR longue/complexe, vu ses limitations), et le filtre multimappeurs a changé (MAPQ vs tag `XT:A:R`).

### 2. Comptages bruts par catégorie (Type 1)

Le fichier qui contient les lectures de Type 1 est *.ssDNA_type1.bed, et son nombre de lignes correspond au nombre de fragments de Type 1. Il est issu de `ITR_id_bwamem.pl`
Logique de classification (dans ITR_id_bwamem.pl) : Pour chaque paire de lectures, le script :

- Analyse le CIGAR et le tag MD de la lecture R2 pour identifier la région candidate de l'ITR (la partie en "soft-clip" et les mismatches en début d'alignement).
- Compare cette région de R2 au complément de la séquence de R1 pour déterminer la longueur du "fill-in" (la partie complémentaire à R1 mais non-homologue au génome).
- Mesure la micro-homologie résiduelle.
- Applique les seuils de classification : une paire est classée Type 1 si elle remplit les critères stricts définis dans la publication originale (et codés dans le script) :
    - ITR_len > 5
    - aln_offset > 2 (longueur du "fill-in")
    - uH_len >= 0 (micro-homologie)

| Pipeline | ssDNA_type1 (T1) |
|---|---|
| Ancien | 6 285 |
| Nouveau | 9 308 |

### 3. Recouvrement génomique (bedtools intersect, `-f 0.5`)

**Type1 ancien -> où il atterrit dans le nouveau pipeline (6 285 lectures) :**

| Catégorie nouveau pipeline | Nombre | % |
|---|---|---|
| Type1 (confirmé) | 5 391 | 85,8 % |
| Type2 | 0 | 0 % |
| dsDNA | 28 | 0,4 % |
| **Unclassified** | **526** | **8,4 %** |
| Absent de toutes catégories (probablement filtré en amont) | 340 | 5,4 % |

**Type1 nouveau -> d'où il vient dans l'ancien pipeline (9 308 lectures) :**

| Catégorie ancien pipeline | Nombre | % |
|---|---|---|
| Type1 (confirmé) | 5 391 | 57,9 % |
| Type2 | 26 | 0,3 % |
| dsDNA | 0 | 0 % |
| Unclassified | 1 | 0 % |
| Absent de toutes catégories (nouvelles lectures récupérées) | 3 890 | 41,8 % |

## Interprétation

- Le gros de l'écart T1 (3 890/3 917 lectures "nouvelles") s'explique par le volume d'entrée, pas par la logique de classification : ces lectures n'existaient dans aucune catégorie de l'ancien pipeline (probablement jamais alignées par `bwa_ra`). Cohérent avec une sensibilité d'alignement améliorée. Résultat plutôt positif
- 526 lectures (8,4 % des T1 de l'ancien pipeline) passent en `unclassified` dans le nouveau 

## Analyse des 526 lectures Type1 → Unclassified

Cette section documente l'investigation complète, y compris une hypothèse initiale qui s'est révélée fausse; gardée ici parce que le raisonnement qui l'a réfutée est ce qui a permis d'arriver à la vraie explication.

### Premier exemple regardé et hypothèse initiale (réfutée)

Exemple tracé : `chr1:3196703` (lecture `E00526:211:HLLMKCCXY:1:1102:22201:12806`).

- **Ancien pipeline** : `ITR_len=15, uH_len=10, aln_offset=5` → Type1.
- **Nouveau pipeline** : R2 alignée `45M5S` par `bwa mem` (0 mismatch sur les 45M) → `ITR_len` calculé = 5 (clip) + 0 (pas de mismatch en tête du MD) = **5**, juste sous le seuil `>5` → `unclassified`.

Hypothèse initiale : la séquence de R2 étant répétitive/pauvre en complexité (motifs `ATAT` récurrents), `bwa mem` aurait trouvé par coïncidence un alignement sans mismatch sur une portion qui serait en réalité de l'ITR, sous-estimant `ITR_len`. Une comparaison directe de R2 à R1 sur 20pb montrait une concordance de 80-100 %, ce qui semblait confirmer que la vraie ITR s'étendait au-delà des 5pb clippées.

une forte concordance à R1 sur une séquence pauvre en complexité ne prouve rien, car ce type de séquence a une probabilité élevée de "bien matcher" n'importe quel comparateur court (génome, R1, autre) simplement par manque d'information. exactement le phénomène qui a fait que `bwa mem` a lui-même trouvé un alignement génomique à cet endroit. La concordance à R1 ne distingue donc pas "vraiment de l'ITR" de "coïncidence de séquence pauvre en complexité"

### Test de l'hypothèse "séquence répétitive" sur l'ensemble des 526 (réfuté)

Comparaison de l'entropie de Shannon (composition en bases, 20 premières bases de R2) entre le groupe désaccord (526 lectures) et un groupe contrôle (lectures Type1 confirmées dans les deux pipelines, même taille d'échantillon) :

| Groupe | Entropie moyenne | n |
|---|---|---|
| Désaccord (T1 ancien → unclassified nouveau) | 1,804 | 526 |
| Contrôle (T1 confirmé dans les deux) | 1,781 | 669 |

Aucune différence significative, et la distribution du groupe désaccord est majoritairement en haute entropie (1,6-2,0), avec seulement une petite traîne de cas répétitifs. L'hypothèse "séquence répétitive" est donc écartée 

### Le vrai motif : CIGAR/MD sur un échantillon plus large

Extraction de CIGAR + MD pour 20 lectures R2 du groupe désaccord, dans les deux pipelines :

- **Nouveau pipeline (`bwa mem`)** : systématiquement `CIGAR` du type `NM/(50-N)S` avec `N` entre 45 et 47 (donc clip de 3 à 5pb), et **zéro mismatch** dans la portion alignée (`MD` = nombre pur, ex. `MD:Z:45`).
- **Ancien pipeline (`bwa_ra`)** : systématiquement `CIGAR=50M` (aucun clip, alignement forcé sur toute la longueur), avec seulement **2 à 5 mismatches explicites** près de la même frontière (ex. `MD:Z:45G0A0A0T1`).

Les deux aligneurs, indépendamment, s'accordent sur seulement 3-5pb de désaccord réel avec le génome à cette frontière jamais 9 ni 15pb. Comme `ITR_len=3`, `4` ou `5` tombe exactement dans la zone grise intentionnelle du schéma de classification (`dsDNA` exige `<3`, `T1/T2` exige `>5`), ces lectures atterrissent en `unclassified` par construction, pas par bug de seuil.

### Reproduction fidèle de `find_offset()` (ancien pipeline) ,cause racine confirmée

Pour trancher définitivement, `find_offset()` a été ré-implémenté fidèlement en Perl et exécuté sur les données réelles de l'exemple `chr1:3196703` (R1 : `50M`, `MD:Z:50`, 0 mismatch ; R2 : `50M`, `MD:Z:45G0A0A0T1`, 4 mismatches). Le résultat reproduit exactement `ITR_len=15, aln_offset=5, uH_len=10` , validant la fidélité de la ré-implémentation.

L'inspection des valeurs intermédiaires révèle le mécanisme exact :

- **`cmpAB`** (R2 comparé au **génome**, reconstruit depuis le tag MD) : seulement **4 positions** sur 50 sont de vrais désaccords avec le génome (`-2`, correspondant exactement aux 4 mismatches du MD). Toutes les autres positions concordent avec le génome (`1`).
- **`seq_offset`** (devenu `ITR_len=15` chez l'appelant) provient d'un calcul séparé et indépendant : R2 comparé au proxy R1 (`arrC`), sans jamais vérifier si les positions concordantes avec R1 concordent *aussi* avec le génome. Ce calcul trouve 15 positions de concordance avant 2 mismatches , mais 11 de ces 15 positions sont des positions où R2 correspond aussi parfaitement au génome (`cmpAB=1`). Elles ne gonflent `ITR_len` que parce que R1 et le génome partagent, à cet endroit, une composition répétitive similaire.

**Conclusion : l'ancien algorithme a un vrai biais méthodologique** : son calcul d'`ITR_len` (`seq_offset`) ne croise jamais la concordance à R1 avec la concordance (ou non) au génome, et peut donc gonfler l'estimation d'ITR sur des régions où R1 et le génome se ressemblent localement, indépendamment de toute vraie ITR. C'est exactement le problème qu'on soupçonnait initialement du côté du nouveau pipeline, mais qui s'applique en réalité à l'ancien. Le nouveau pipeline (`ITR_len=5`, borné par ce que `bwa mem` a jugé non-génomique) est méthodologiquement plus défendable sur ce cas précis.

## Statut / suite

- **Conclusion révisée** : les 526 lectures Type1→Unclassified sont vraisemblablement d'anciens faux positifs de `find_offset()` (ancien pipeline) corrigés par la nouvelle approche, plutôt que des régressions du nouveau pipeline. Le nouveau `find_offset_mem()` n'a pas montré de bug dans cette investigation. au contraire, sa logique (borner la région ITR candidate par ce qu'un aligneur pangénomique a jugé non-génomique, puis seulement ensuite comparer à R1) s'est révélée plus robuste que celle de l'ancien algorithme face aux régions localement répétitives.
- Piste d'amélioration envisagée puis écartée : étendre la recherche de fill-in au-delà de la frontière de clip de `bwa mem` en utilisant l'identité à R1 comme signal principal. Abandonnée après investigation : cette approche réintroduirait précisément le biais identifié dans l'ancien algorithme (confondre "concorde avec R1" et "est de l'ITR" sans vérifier la concordance génomique).
- Impact estimé : ~8 % des Type1 historiques déclassés en `unclassified` par le nouveau pipeline, concentré sur les frontières ITR courtes (3-5pb, zone grise intentionnelle du schéma de classification) plutôt que sur des séquences répétitives comme on le pensait initialement.
- Reste à faire (partiellement fait, voir section suivante) : même comparaison sur données NovaSeq réelles une fois qu'un jeu de données NovaSeq *avec ancien pipeline fonctionnel* sera disponible (impossible actuellement puisque l'ancien pipeline plante justement sur NovaSeq).

## Run NovaSeq réel (post-migration, sans comparaison possible)

**Date** : 2026-09-15. **Contexte** : premier run de `main.nf` sur données NovaSeq réelles de production (échantillon `P9_HET`, 2 réplicats IP `R1`/`R2` + 1 contrôle input `IN_R1`, ~330-340M reads chacun) pas le petit jeu de test. Terminé de bout en bout sans erreur : le crash `bwa sampe` d'origine est bien résolu.

Comme il est structurellement impossible d'avoir une référence ancien-pipeline sur NovaSeq (c'est exactement ce qui plantait), pas de comparaison ancien/nouveau possible ici. La vérification s'est donc limitée à un contrôle de plausibilité des sorties (`flagstat`, `SSDSreport.tab`), pas une validation de fidélité à l'algorithme d'origine comme pour le HiSeq ci-dessus.

### Flagstat (bwaAlign)

| Échantillon | Mapped | Properly paired | Duplicats |
|---|---|---|---|
| R1 (IP) | 92,37 % | 66,19 % | 4,35 % |
| R2 (IP) | 93,85 % | 65,22 % | 4,25 % |
| IN_R1 (input) | 95,41 % | 67,74 % | 3,65 % |

Le taux de "properly paired" (~65-68 %), plus bas qu'un ChIP-seq classique, est identique entre les 2 réplicats IP **et** le contrôle input, cohérent avec une caractéristique du run/de la prep de librairie plutôt qu'un artefact d'alignement spécifique à l'ITR/fill-in.

### Répartition des catégories SSDS (`SSDSreport.tab`, `totinfo`)

| Catégorie | R1 | R2 |
|---|---|---|
| Type1 | 65,4 % | 68,8 % |
| Type2 | 0,16 % | 0,13 % |
| dsDNA | 14,1 % | 10,3 % |
| unclassified | 20,2 % | 20,7 % |

Majoritairement Type1 dans les deux réplicats : signature attendue d'un ChIP-DMC1 qui fonctionne. Les histogrammes `ITR_ssDNA_type1`/`Fragments_ssDNA_type1` décroissent de façon lisse (pas de pic dégénéré), cohérent avec une vraie distribution physique de fill-in plutôt qu'un artefact de parsing.

FRIP (vs pics de novo) : R1 (T1=0,214 ; T2=0,202 ; dsDNA=0,098/0,072 ; unclassified=0,230) systématiquement plus bas que R2 (T1=0,292 ; T2=0,273 ; dsDNA=0,166/0,128 ; unclassified=0,297), non pathologique dans les deux cas (ni 0 ni 1), mais voir l'écart ci-dessous.

### Écart de fragments uniques Type1 entre réplicats élucidé via `satCurve.tab`

R1 et R2 ont un nombre brut de fragments Type1 comparable (63,3M vs 65,2M) mais un nombre de fragments **uniques** Type1 très différent : 23,7M (37 % d'uniques) pour R1 contre seulement 8,55M (13 % d'uniques) pour R2, alors que le taux de duplicats au niveau BAM global est quasi identique entre les deux (4,35 % vs 4,25 %). L'écart est donc spécifique à la catégorie Type1, pas une différence générale de duplication.

`02_results/peaks/with-input/saturation_curve/standard/satCurve.tab` (nombre de hotspots `hs` appelés par MACS2 à différents paliers de sous-échantillonnage des reads Type1) tranche entre les deux hypothèses envisagées :

| Palier | Reads R1 | hs R1 | Reads R2 | hs R2 | ratio hs R2/R1 |
|---|---|---|---|---|---|
| 20 % | 12,6M | 49 485 | 12,9M | 82 349 | 1,66× |
| 40 % | 25,3M | 44 429 | 25,7M | 70 256 | 1,58× |
| 60 % | 37,9M | 38 541 | 38,6M | 80 360 | 2,08× |
| 80 % | 50,5M | 41 833 | 51,5M | 89 078 | 2,13× |
| 100 % | 63,2M | 40 519 | 64,3M | 79 401 | 1,96× |

À chaque palier de reads comparable, R2 appelle systématiquement ~1,6-2,1× plus de hotspots que R1, cohérent avec l'écart de FRIP déjà observé (R2 nettement au-dessus de R1 sur toutes les catégories). **Conclusion** : R2 est une librairie moins complexe (moins de molécules uniques d'origine, d'où le taux de fragments uniques Type1 plus faible) mais plus concentrée/mieux enrichie sur les vrais hotspots (d'où le FRIP et le nombre de hotspots plus élevés à profondeur égale). une différence de qualité de librairie entre réplicats au niveau paillasse, pas un artefact du pipeline ni de la migration bwa mem. Point de QC à garder en tête pour l'interprétation biologique (comparaison des réplicats), mais ne remet pas en cause la migration.



### Hypothèse chimie 2-couleurs NovaSeq (G-tailing) testée sur les données disponibles

Hypothèse : la chimie 2 couleurs de NovaSeq (un cycle sans signal peut être appelé "G" au lieu de "N", contrairement à HiSeq en 4 couleurs) pourrait produire un artefact spécifiquement en bout de read  potentiellement pile là où se joue la classification ITR (5' de R2). Vérifiée via les données FastQC/MultiQC déjà synchronisées pour le contrôle input `IN_R1` (`02_results/qc/multiqc/.../multiqc_fastqc.txt`, `mqc_fastqc_per_sequence_gc_content_plot_Percentages.txt`, `mqc_fastqc_adapter_content_plot_1.txt`) :

- **Read1** (biologiquement R1, censé être 100 % génomique) : ~20 % des séquences tombent dans un pic anormal à 92-99 % GC (le pic génomique normal culmine ~40 % et retombe à ~0 dès 55-60 % GC) ; `per_base_sequence_content` = **fail** ; séquence la plus fréquente = 5,77 % de toutes les reads. Fichier vérifié = **post-trimming complet** (`crop_trim_R1.fastq.gz`). l'artefact traverse donc le trimming actuel sans être nettoyé. Le module "adapter content" de FastQC ne remonte aucune courbe pour Read1 (sous seuil) : pas de l'adaptateur classique, cohérent avec un vrai artefact de dark-cycle plutôt qu'un read-through d'adaptateur.
- **Read2** (biologiquement R2, celui qui porte l'ITR/fill-in) : même signature mais bien plus faible (~0,2 % de densité contre ~20 % sur Read1) ; en revanche Read2 montre une vraie courbe d'adaptateur universel croissante avec la position. cohérent avec le phénomène "adapter" déjà expliqué plus haut (lecture à travers un insert court), un phénomène distinct du G-tailing.

**Conclusion** : le G-tailing NovaSeq est réel et mesurable dans ce jeu de données (pas une fausse piste), mais concentré sur Read1, pas sur Read2 où se joue la classification ITR. risque faible pour la fidélité de la classification. Une read quasi tout-G ne s'aligne de toute façon pas correctement : elle sort vraisemblablement dans la fraction non-mappée du flagstat (~5-8 %) plutôt que de fausser silencieusement une classification.

**Point actionnable indépendant trouvé au passage** : le process `trimming` de `main.nf` appelle `trim_galore` sans l'option `--nextseq` (mode dédié à la chimie 2 couleurs NovaSeq/NextSeq, qui traite les runs de G terminaux comme du bruit à couper). Amélioration naturelle à coupler avec `--sequencer novaseq`, sur le même principe que l'ajustement déjà fait pour `OPTICAL_DUPLICATE_PIXEL_DISTANCE`. Pas encore implémenté.

**Réserve** : vérifié uniquement sur le contrôle input (FastQC des réplicats IP R1/R2 pas encore synchronisé dans ce rapport multiqc au moment du test) ; analyse sur stats agrégées (%GC, courbe adaptateur), pas sur le texte des séquences lui-même. "cohérent avec" un artefact de dark-cycle, pas confirmé au niveau séquence.

### Conclusion à ce stade

Rien de dégénéré ou d'alarmant sur les métriques inspectées. le run NovaSeq est plausible et exploitable. Mais "tourne proprement" et "distribution plausible" ne valent pas fidélité algorithmique démontrée (contrairement au HiSeq, impossible à établir ici faute de référence). L'écart Type1 unique R1/R2 est élucidé (différence de complexité/enrichissement de librairie entre réplicats, pas un problème de pipeline)
Reste ouvert : l'inspection du contrôle input au-delà du flagstat 
