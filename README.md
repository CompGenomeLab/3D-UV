
- XR-seq and Damage-seq analysis

Raw XR-seq and Damage-seq datasets are processed using pipeline in [here](https://github.com/CompGenomeLab/xr-ds-seq-snakemake)

- TAD calling


TAD calling is performed using Arrowhead with default parameters using script in `tad-calling.sh`

- 3D modelling


TADs and contact matrices are first processed with `3d-preprocess.sh` by following procedure in [here](https://www.nature.com/articles/nprot.2018.009).
Then Chrom3D used for generating 3D genome models via `3d-model.sh`

- Downstream analysis


Codes in `uv-3d.py` is used for analysis of Damage-seq and XR-seq data on produced 3D genome model.
