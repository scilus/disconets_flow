# Disconects_flow

Disconects_flow allows you to know the impact of a cavity on structural connectivity matrices.

Yeo atlas used in the paper [link](https://box.criugm.qc.ca/f/65e07378c3374453ae9c/?dl=1)

Subset of tractograms [link](https://box.criugm.qc.ca/f/034b74d2c9844da38951/?dl=1)

### Command line:
```
nextflow run main.nf \
        --input [FullPathTo]data/ \
        --atlas [FullPathTo]/atlas/ \
        --tractograms [FullPathTo]/tractograms/ \
        -with-singularity scilus-1.0.0-rc1.sif \
        resume
```


### Disconects input structure
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
