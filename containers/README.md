# Singularity Containers

## python-3.8.sif
```
singularity pull python-3.8.sif docker://python:3.8-bullseye
```

Source : https://hub.docker.com/layers/library/python/3.8-bullseye/images/sha256-c6102cb482b56fe46114f6087c8286d0d2470c16cc16e8a060164a95dac19ac5?context=explore   
   

## bioinformatics.sif
```
singularity pull bioinformatics.sif shub://thakk/bioinformatics:latest
```

Source : https://singularityhub.github.io/singularityhub-archive/containers/thakk-bioinformatics-latest/

## trimming-box-1.0.sif
# cf https://hub.docker.com/layers/pauffret/trimming-box/1.0/images/sha256-c24171d2f9869c0783b966e4dbb6a143a15543bb6c8fa790c6923d820f730ee1?context=repo
# In trimming-box directory containing only Dockerfile and environment.yml :
```
sudo docker build -t pauffret/trimming-box:1.0 .

sudo docker push pauffret/trimming-box:1.0

singularity pull trimming-box_1.0.sif docker://pauffret/trimming-box:1.0
```

## bam-box-1.0.sif
# cf https://hub.docker.com/layers/pauffret/bam-box/1.0/images/sha256-d356d4c07d35feea8d39209f94c6047b4cbd6bb1a306b477835191affe1de83a?context=repo
# In bam-box directory containing only Dockerfile and environment.yml :
```
sudo docker build -t pauffret/bam-box:1.0 .

sudo docker push pauffret/bam-box:1.0

singularity pull bam-box_1.0.sif docker://pauffret/bam-box:1.0
```

