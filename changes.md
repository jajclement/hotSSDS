modifications apportées : 

# Migration bwa aln/bwa-ra/sampe → bwa mem, et ajustements NovaSeq

## Contexte

Sur données NovaSeq, l'étape `bwa sampe` du pipeline plantait de façon reproductible :

```
[bwa_sai2sam_pe_core] refine gapped alignments...
bwa_0.7.12: bwase.c:180: bwa_refine_gapped_core: Assertion `re - rb == rlen' failed.
```

Il s'agit d'un bug connu de `bwa` 0.6/0.7 (`bwa_refine_gapped_core`), vraisemblablement aggravé par le binning des scores de qualité NovaSeq (RTA3, ~4 valeurs de qualité au lieu d'un gradient continu), qui multiplie les cas d'égalité de score lors du raffinement gapped.

Le pipeline utilisait un alignement en trois étapes propre à ce projet :
- R1 aligné avec `bwa aln` (binaire stock `bin/bwa_0.7.12`)
- R2 aligné avec `bwa_ra aln` (binaire modifié `bin/bwa_ra_0.7.12`, dit "bwa right-align", conçu par K. Brick en 2012 pour forcer l'alignement à travers la partie non-homologue en 5' de R2, l'ITR)
- fusion des deux avec `bwa sampe` (binaire stock) c'est cette étape qui plante

Aucune source du binaire `bwa_ra` n'a pu être retrouvée publiquement (vérifié dans `kevbrick/SSDSpipeline`, `kevbrick/SSDSnextflowPipeline`, et le profil GitHub complet de l'auteur), il n'était donc pas possible de corriger ou recompiler ce binaire.

Le papier de méthode (Khil et al. 2012, *Genome Research*, section "Computational ssDNA detection in the sequencing libraries") précise que le rôle de `bwa-ra` est uniquement de positionner R2 sur le génome malgré l'ITR, et autorise explicitement un autre aligneur : *"any other aligner can be used provided that they will find the longest 3′ end subsequence for the second read."* `bwa mem` fait exactement cela nativement via le soft-clip, sans bug de raffinement gapped puisqu'il ne repose pas sur l'étape `sampe`.




## 1. Alignement : `bwa aln` + `bwa-ra aln` + `bwa sampe` remplacés par un seul appel à `bwa mem`

**Fichier : [main.nf](main.nf), process `bwaAlign`**

- Les trois commandes d'alignement sont remplacées par un seul appel :
  ```
  ./${bwa} mem -M -t ${task.cpus} ${fasta} ${fqR1} ${fqR2} > ${tmp_name}.unsorted.sam
  ```
  Aucun nouveau binaire n'a été nécessaire : `bin/bwa_0.7.12` supporte déjà la sous-commande `mem` (introduite dans bwa 0.7.0). `bwa mem` gère le pairing nativement, donc l'étape `sampe` et son bug disparaissent.
- Le filtre multimappeurs (`--with_multimap false`) reposait sur le tag `XT:A:R`, spécifique à `bwa aln` et absent de `bwa mem`. Il a été remplacé par un filtre MAPQ (`samtools view -q 1`, `bwa mem` rapportant MAPQ 0 pour les reads multi-mappés).
- Le process ne prend plus `bwa_ra` en entrée ; les channels `bwa_ra_ch` et `bwa_all_ch` ont été supprimés.
- références orphelines `bwa_ra` / `--custom_bwa_ra` 
- Le binaire `bin/bwa_ra_0.7.12` a été laissé sur disque 




## 2. Nouveau script de classification ITR : `bin/ITR_id_bwamem.pl`

**Nouveau fichier : [bin/ITR_id_bwamem.pl](bin/ITR_id_bwamem.pl)**, remplace `bin/ITR_id_v2c_NextFlow2.pl` dans le process `parseITRs` (`main.nf`).

### Nouveau script ITR_id_bwamem.pl

L'ancien script (`find_offset()`) supposait la structure d'alignement produite par `bwa aln`/`bwa-ra`/`sampe`, où l'ITR de R2 était forcée dans l'alignement comme des mismatches positionnés, jamais comme un vrai soft-clip. Avec `bwa mem`, l'ITR est représenté :
- en <b>soft-clip** (CIGAR `S`)</b> pour une ITR longue 
- en <b>mismatches résiduels** dans le bloc `M` (via le tag `MD`)</b> pour une ITR courte 

L'ancien algorithme n'aurait pas su interpréter cette nouvelle structure sans produire une classification T1/T2/dsDNA fausse

### Principe du nouvel algorithme (`find_offset_mem`)

Reprend la définition du papier (*"The fill-in part of the ITR must be complementary to the 5'-end of the first read but not homologous to the target genomic DNA"*) :

1. Lit le CIGAR et le tag `MD` de R2 pour extraire la région ITR candidate : longueur du soft-clip 5' **+** mismatches consécutifs en tête du bloc `M` (gère l'orientation du brin, le soft-clip du 5' biologique pouvant être stocké à droite du CIGAR/`SEQ` si R2 est sur le brin `-`).
2. Compare cette région, base par base, à la séquence de R1 elle-même (identité directe, pas complément), en tolérant 1 mismatch donne la longueur du "fill-in".
3. Étend la comparaison au-delà pour mesurer la microhomologie 
4. Retourne `(ITR_len, aln_offset, uH_len, aln_correction)` avec la même signature que l'ancien `find_offset()`.

Tout le reste du script (ouverture des fichiers BED/SAM, seuils de classification T1/T2/dsDNA/dsDNA_strict, réécriture du CIGAR/MD/NM via `generateSAMline`/`trimCigar`/`trimXM`) est repris à l'identique, la sortie (noms de fichiers) est inchangée, donc `parseITRs` n'a pas eu besoin d'être modifié au-delà du nom du script appelé.

### Limite connue, à valider

Une insertion CIGAR tombant exactement à la frontière ITR/génome n'est pas prise en compte lors de la lecture du tag `MD` (le `MD` n'encode pas les insertions). Ultra rare et non bloquant.

### Bugs trouvés et corrigés sur données réelles

Deux bugs bloquants ont été découverts en rejouant le pipeline sur un jeu de test réel, invisibles sur les tests unitaires initiaux car ceux-ci partageaient la même hypothèse erronée que le code testé :

1. **Sens de comparaison inversé dans `find_offset_mem()`** : la version initiale comparait R2 au **complément** de R1 au lieu de R1 directement (erreur de sens signalée ci-dessus, point 2). Conséquence : `aln_offset` plafonnait quasi systématiquement à 0-2, donc aucune paire n'était jamais classée en Type1. Diagnostiqué en rejouant la fonction sur une vraie paire de reads et en observant qu'ajouter un complément supplémentaire donnait une concordance parfaite sur toute la longueur de l'ITR : signal sans ambiguïté. Corrigé en comparant R2 directement à R1 (identité).
2. **Coordonnées BED invalides pour les fragments courts** : `$aln_correction` pouvait dépasser la moitié de l'étendue du fragment, inversant `start`/`end` et faisant planter `genomeCoverageBed`. Corrigé en plafonnant `$aln_correction` à la moitié de l'étendue du fragment.

Voir [VALIDATION_bwa_mem_vs_original.md](VALIDATION_bwa_mem_vs_original.md) pour la comparaison complète ancien vs nouveau pipeline, qui a aussi révélé que l'ancien algorithme `find_offset()` a son propre biais (surestimation de l'ITR sur les régions où R1 et le génome se ressemblent localement).




## 3. Paramètre `--sequencer` (hiseq/novaseq) pour la distance de duplicats optiques

**Fichier : `main.nf`** : ajouté indépendamment de la migration `bwa mem`, pour un problème NovaSeq différent (sans lien avec le bug `sampe`).

- Nouveau paramètre `--sequencer` (`'hiseq'` par défaut, ou `'novaseq'`) 
- Le paramètre Picard `OPTICAL_DUPLICATE_PIXEL_DISTANCE` (utilisé dans les 3 appels `MarkDuplicatesWithMateCigar` de `bwaAlign` et `filterBam`) est calculé dynamiquement :
  - si `--sequencer novaseq` et `--picard_optdup_distance` est resté à sa valeur par défaut (100, calibrée pour les flowcells non-patternés type HiSeq2500), il est automatiquement ajusté à 2500 (valeur recommandée pour les flowcells patternés type NovaSeq) ;
  - si l'utilisateur a explicitement fourni une autre valeur avec `--sequencer novaseq`, un avertissement est émis mais sa valeur est respectée ;
  - avertissement également si `--sequencer` a une valeur inconnue (ni `hiseq` ni `novaseq`).
- Documenté dans l'aide CLI de `main.nf` et dans le README.

## 4. `--sequencer novaseq` bascule aussi `trim_galore` en mode `--nextseq`

**Fichier : `main.nf`, process `trimming`** : trouvé en creusant un run NovaSeq réel de production (voir [VALIDATION_bwa_mem_vs_original.md](VALIDATION_bwa_mem_vs_original.md), section "Hypothèse chimie 2-couleurs NovaSeq").

la chimie 2 couleurs de NovaSeq peut appeler un cycle sans signal "G" au lieu de "N" (contrairement à HiSeq, 4 couleurs). Confirmé sur données réelles : ~20 % des lectures R1 d'un échantillon tombaient dans un pic anormal à 92-99 % GC, non nettoyé par le trimming existant (mesuré **après** trim_galore) et non détecté comme de l'adaptateur classique par FastQC.
- `trim_galore` dispose d'un mode dédié (`--nextseq <seuil>`) qui ignore la qualité des runs de G terminaux au lieu de les traiter comme des bases de bonne qualité. Il remplace `--quality` (les deux ne sont pas censés être combinés).
- Comportement : si `--sequencer novaseq`, le process `trimming` utilise `--nextseq ${params.trimg_quality}` à la place de `--quality ${params.trimg_quality}` dans les deux branches `trim_galore` (avec ou sans fichier d'adaptateurs fourni). Sinon (`hiseq`, valeur par défaut), comportement inchangé.
- Ne s'applique qu'avec `--with_trimgalore` : `trimmomatic` n'a pas d'équivalent au mode `--nextseq`. Un avertissement est émis si `--sequencer novaseq` est utilisé sans `--with_trimgalore`.
- Risque : l'analyse sur le run NovaSeq a montré que l'artefact touche surtout R1 et qu'une lecture quasi tout-G ne s'aligne simplement pas (sort dans les reads non mappées plutôt que de fausser une classification). Correctif préventif de qualité de données, pas un correctif de bug de classification.

