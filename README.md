# Disconets_flow

Disconets_flow allows you to analyze the impact of a cavity on structural connectivity matrices.

Yeo atlas used in the paper [link](https://box.criugm.qc.ca/f/65e07378c3374453ae9c/?dl=1).

Please cite:
```
Mrah S, Descoteaux M, Wager M, Boré A, Rheault F, Thirion B, Mandonnet E. Network-level prediction of set-shifting deterioration after lower-grade glioma resection. J Neurosurg. 2022 Mar 4:1-9. doi: 10.3171/2022.1.JNS212257. Epub ahead of print. PMID: 35245898.
```

### Build singularity or docker image
```
# Singularity
singularity build scilus_2.0.2.sif docker://scilus/scilus:2.0.2

# Docker
docker pull scilus/scilus:2.0.2
```

## Run Disconets_flow
```
# With singularity image
nextflow run main.nf \
        --input [FullPathToLesions]/ \
        --atlas [FullPathToAtlas]/ \
        --tractograms [FullPathToTractograms]/ \
        -with-singularity scilus_2.0.2.sif \
        resume

# With docker image
nextflow run main.nf \
        --input [FullPathToLesions]/ \
        --atlas [FullPathToAtlas]/ \
        --tractograms [FullPathToTractograms]/ \
        -with-docker scilus/scilus:latest \
        resume
```


### Disconets input structure
```
Root
├── S01
│   ├── *t1.nii.gz (optional)
│   └── *cavity.nii.gz
└── S02
    ├── *t1.nii.gz (optional)
    └── *cavity.nii.gz
 ```

 ### Tractograms input structure
 ```
 Root
 ├── T01
 │   ├── *t1.nii.gz (optional)
 │   └── *.trk
 └── T02
     ├── *t1.nii.gz (optional)
     └── *.trk
  ```
