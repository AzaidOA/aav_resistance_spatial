# AAV Resistance Signature Mapping in Mouse Brain

This repository contains the computational pipeline and associated scripts to identify and spatially map microenvironmental barriers to AAV-mediated gene delivery in mouse brain tissue, specifically in the context of prion disease gene therapy.

## Project Overview

The goal of this project is to understand why AAV transduction efficiency varies across neurons in the cortex and hippocampus, using a model inspired by the zinc finger-based PRNP knockdown therapy. By integrating single-nucleus RNA-seq data (scAAVengr) with spatial transcriptomics, we aim to:

- Define gene expression signatures of AAV-resistant cells.
- Project those signatures onto spatial brain data (10x Visium or equivalent).
- Identify anatomical and microenvironmental features (e.g., ECM, immune activation) that correlate with low transduction.
- Make actionable suggestions for improving gene therapy delivery (e.g., enzyme co-delivery).

## Data Sources

- **scRNA-seq (AAV9 proxy):**  
  Öztürk et al., *eLife* 2021 — GSE155242  
  Cortex and hippocampus, targeted and non-targeted AAV enrichment.

- **Spatial transcriptomics:**  
  10x Genomics public mouse brain datasets (cortex and hippocampus).

