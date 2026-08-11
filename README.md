# Reproducibility materials for a cross-disease endometrial omics framework

**A cross-disease endometrial omics framework maps heterogeneous timing and context-specific molecular remodeling in implantation-associated disorders**

## Repository purpose

This repository contains selected public source tables, analysis scripts, metadata, environment summaries and selected source-table-derived outputs supporting the reported analyses. Release v1.1.0 extends the previously archived materials with independent external validation of the frozen analytical timing coordinate.

## Contents

- `source_data/`: selected verified source tables and outputs;
- `source_data/external_validation/`: row-level timing scores, Supplementary Table S15 source data, frozen-gene coverage and validation provenance for four independent external cohorts;
- `scripts/`: selected historical and verified rebuild scripts;
- `scripts/external_validation/`: lightweight reconstruction scripts for the released external-validation score tables;
- `metadata/`: numbering crosswalks, public accessions and formal citations;
- `figure_source_data_index/`: figure-to-source mapping;
- `additional_file_3_moved/`: selected supplementary-table source files retained from the earlier reproducibility package;
- `environments/`: runtime and package summaries;
- `outputs/`: selected source-table-derived figure outputs;
- `documentation/`: reproducibility documentation, including the external-validation scope and interpretation boundaries.

## Public data accessions

- `GSE250130`
- `GSE111974`
- `GSE58144`
- `GSE287278`
- `GSE179640`
- `GSE213216`
- `GSE214411`
- `GSE135485`
- `GSE263897`
- `GSE244236`
- `GSE190580`
- `GSE157718`
- `GSE78851`
- `E-MTAB-14039`
- `GSE234368`
- `GSE98386`
- `GSE180485`
- `GSE144895`

These are public-data references only; raw data are not mirrored here. Original public repositories remain the authoritative source for raw expression data.

## Independent external validation added in v1.1.0

The frozen analytical timing coordinate was evaluated in four additional public cohorts without gene reselection, model refitting or target-cohort recalibration. The released materials support Figure 7, Supplementary Figures S22–S25 and Supplementary Table S15.

Detailed scope, cohort-level inference units and interpretation boundaries are documented in `documentation/EXTERNAL_VALIDATION_README.md`.

## Scope limitation

This repository is a selected source-table reproducibility package supporting the reported analyses.

It is **not**:

- a mirror of the raw GEO datasets;
- a complete collection of large expression matrices;
- a complete archive of all single-cell integration objects;
- a full end-to-end raw-GEO-to-atlas reconstruction pipeline;
- a clinically deployable timing, diagnostic or prediction tool.

Basic inspection should start with `metadata/FINAL_PUBLIC_NUMBERING_CROSSWALK.tsv`, `figure_source_data_index/`, `source_data/external_validation/`, source tables, reports and environment summaries. Re-running selected scripts may require separately obtaining the public accessions and matching the documented environment.

## Citation and contact

Please cite the associated manuscript and the relevant archived release when using these materials. Repository release version: `1.1.0`. GitHub repository: https://github.com/yueminhou22-source/implantation-failure-atlas.

## Archival DOI

The earlier v1.0.0 reproducibility release is archived on Zenodo.

Version-specific DOI (v1.0.0):
https://doi.org/10.5281/zenodo.21849398

Concept DOI (all versions):
https://doi.org/10.5281/zenodo.21849399

The version-specific DOI for v1.1.0 will be added here after the new Zenodo version is published.

Authors: Qing Gao; Tingting Zhang; Wei Liu; Tiantian Ji; Dan Zou; Yuemin Hou.

Department of Obstetrics and Gynecology, The Second Affiliated Hospital of Xi'an Jiaotong University, Xi'an, Shaanxi 710004, China.

Corresponding author: Yuemin Hou — `18896500386@163.com`.

## License

Software code and scripts in this repository are licensed under the MIT License.

Research source tables, metadata, documentation, reports, and other non-software research materials are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

See `LICENSE.md` and the files in `LICENSES/` for details.
