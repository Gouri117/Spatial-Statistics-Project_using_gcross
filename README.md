# Spatial Interaction Analysis of the Colorectal Cancer Tumor Microenvironment

## Overview

This project applies G-cross spatial statistics to characterize cell–cell interaction patterns in colorectal cancer (CRC) tissue. Using spatially resolved single-cell data, we quantify directional spatial relationships between cell types and evaluate how these interaction patterns associate with clinical outcomes.

The focus is on capturing spatial organization and interaction structure within the tumor microenvironment rather than cell-type abundance alone.

Data Description

Data type: Spatial imaging–derived single-cell maps (e.g., multiplex IHC / CODEX)

Resolution: Single-cell coordinates

Annotations include:

Cell type labels

- Spatial (x, y) coordinates

- Sample or image identifiers

- Clinical outcome metadata

Raw imaging data are not included due to data-sharing restrictions.

### Analysis Workflow
1. Spatial Preprocessing

Parsing cell-level spatial coordinates

Filtering low-quality cells and rare cell types

Aligning spatial data with clinical metadata

2. G-Cross Spatial Statistics

Computation of directional G-cross functions between cell-type pairs

Quantification of spatial attraction and avoidance

Calculation of AUC-based interaction scores per sample

3. Feature Engineering

Construction of spatial interaction feature matrices

Aggregation of interaction scores across tissue regions

Normalization across samples

4. Outcome Association

Statistical association of spatial features with clinical outcomes

Identification of interaction signatures linked to prognosis

Exploratory modeling for outcome prediction

### Key Questions

1. Which cell–cell spatial interactions characterize the CRC tumor microenvironment?

2. How do immune–tumor and stromal–immune interactions vary by outcome?

3. Do spatial interaction features provide information beyond cell-type abundance?

4. Which interaction patterns suggest immune exclusion or immune engagement?
