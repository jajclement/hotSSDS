modifications apportées à `main_dsl2.nf` :

# Portage DSL2 : migration bwa mem/NovaSeq et corrections de bugs

Ce fichier documente les changements apportés spécifiquement à [main_dsl2.nf](main_dsl2.nf). Pour le contexte complet de la migration `bwa mem`/NovaSeq elle-même faite d'abord sur `main.nf` : voir [changes.md](changes.md). 

## 1. Portage de la migration bwa mem/NovaSeq vers `main_dsl2.nf`

Les quatre changements de changes.md (migration `bwa mem`, script `ITR_id_bwamem.pl`, paramètre `--sequencer`, mode `--nextseq`) ont été reportés à l'identique dans `main_dsl2.nf` (process `bwaAlign`, `parseITRs`, `trimming`, et bloc `workflow{}` pour la validation/l'ajustement des paramètres). Portage syntaxique direct, aucune divergence de logique volontaire avec `main.nf` sur ce périmètre.


## 2. Correction de bug bloquants de `main_dsl2.nf`

Voir [DSL2_PORT_REVIEW.md](DSL2_PORT_REVIEW.md)