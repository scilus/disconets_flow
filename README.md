# Disconets_flow

Disconets_flow allows you to analyse the impact of a cavity on structural connectivity matrices.

Yeo atlas used in the paper [link](https://box.criugm.qc.ca/f/65e07378c3374453ae9c/?dl=1)

Subset of tractograms [link](https://box.criugm.qc.ca/f/034b74d2c9844da38951/?dl=1)



### Build singularity or docker image
```
# Singularity
sudo singularity build scilus_latest.sif docker://scilus/scilus:latest

# Docker
sudo docker pull scilus/scilus:latest
```

## Run Disconets_flow
```
# With singularity image
nextflow run main.nf \
        --input [FullPathTo]data/ \
        --atlas [FullPathTo]/atlas/ \
        --tractograms [FullPathTo]/tractograms/ \
        -with-singularity scilus_latest.sif \
        resume

# With docker image
nextflow run main.nf \
        --input [FullPathTo]data/ \
        --atlas [FullPathTo]/atlas/ \
        --tractograms [FullPathTo]/tractograms/ \
        -with-docker scilus/scilus:latest \
        resume
```


### Disconets input structure
```
Input
├── sub-01
│   ├── cavity.nii.gz
└── sub-02
    └── cavity.nii.gz
 ```

 Please cite:
 ```
 Mrah S, Descoteaux M, Wager M, Boré A, Rheault F, Thirion B, Mandonnet E. Network-level prediction of set-shifting deterioration after lower-grade glioma resection. J Neurosurg. 2022 Mar 4:1-9. doi: 10.3171/2022.1.JNS212257. Epub ahead of print. PMID: 35245898.
 ```
