# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

stGEAToolbox is a MATLAB toolbox for **spatial transcriptomics** analysis. It is a companion to scGEAToolbox and **requires scGEAToolbox on the MATLAB path** to function — stGEAToolbox has no standalone mode.

**Entry point:**
```matlab
stgeatool          % Launch the GUI (prompts for data source)
stgeatool(ste)     % Launch directly with a SpatialTranscriptomicsExperiment
```

## Hard Dependency on scGEAToolbox

stGEAToolbox calls scGEAToolbox functions extensively. Key scGEAToolbox functions used:
- `sc_readhdf5file`, `sc_readtsvfile`, `sc_readgeoaccess` — data I/O
- `sc_norm`, `sc_cellscore`, `sc_pickmarkers`, `sc_knngraph` — analysis
- `sc_cluster_s`, `sc_cluster_x` — clustering
- `SingleCellExperiment` class — the inner data container
- `gui.*`, `pkg.*`, `ten.*` packages — GUI utilities, stats, network analysis

**Path resolution pattern for scGEAToolbox assets:**
- `i_ctypescore.m` locates the PanglaoDB marker file via `fileparts(which('scgeatool'))` + `external/fun_alona_panglaodb/marker_hs.mat`
- `st_ridgecci.m` calls `[a, b] = cdgea` to get scGEAToolbox root, then loads `assets/Ligand_Receptor/Ligand_Receptor.mat`

When scGEAToolbox asset paths change, update these two files accordingly.

## Core Data Structure

`SpatialTranscriptomicsExperiment` (`@SpatialTranscriptomicsExperiment/`) wraps a `SingleCellExperiment` with spatial information:

```
ste.sce     → SingleCellExperiment (expression matrix X, gene names g, cell metadata)
ste.xy      → N×2 spatial coordinates of spots
ste.img     → histology image array
ste.tissue_positions_list  → table from 10x spatial/ folder
ste.scalefactors_json      → struct from scalefactors_json.json
ste.metadata               → string array of provenance notes
```

## Architecture

### Package Structure
- `@SpatialTranscriptomicsExperiment/` — Core class (constructor, spatial transforms, accessors)
- `+st/+gui/` — GUI helpers (update checker, cluster selection dialog, menu utilities)
- `+st/+pkg/` — Spatial utilities (e.g. `polycentroid.m` for cell polygon centroids)
- `+st/+run/` — R tool wrappers: `BayesSpace.m`, `STdeconvolve.m`, `SpaTalk.m`
- `+st/+run/external/R_*/` — R scripts (`script.R`, `require.R`) for each wrapped tool
- `+st/+run/private/commoncheck_R.m` — shared R environment validator used by all `+run` wrappers
- `private/countmember.m` — frequency utility

### Root-Level Functions
- **Data loading:** `st_read10xdir`, `st_readsegmenteddir`, `st_readgeoaccession`, `st_readgeoaccessionST`, `st_sce2ste`
- **Analysis:** `st_anova` (spatially variable genes), `st_distgrad` (distance gradients), `st_ridgecci` (cell-cell interactions via ligand-receptor), `st_ridgenet` (network-based CCI via TenifoldXct)
- **Clustering:** `st_cluster_img` (image-based), `st_cluster_spo` (expression-based), `cluster_banksy`, `cluster_resnet18`
- **Visualization:** `st_scatter_ste` (main GUI), `st_showccires` / `st_showccires2` / `st_showccires3` (CCI results)
- **Utilities:** `i_ctypescore` (PanglaoDB cell typing), `i_knndist` (KNN distances), `st_sampleimg`, `st_immidist`

### Main GUI: `st_scatter_ste.m`
This ~19,500-line file is the entire interactive analysis environment. It presents three toolbars:
1. **Spatial Transform** — zoom, rotate, flip, move
2. **Analysis & Display** — dot size/alpha, colormaps, 2D/3D toggle, image overlay, clustering, cell scoring, marker genes, export
3. **Statistics** — SVG detection, multiview linked brushing, PowerPoint export, link to scgeatool

State is stored via `guidata()`. Callbacks are nested functions within `st_scatter_ste`.

### R Tool Integration Pattern
All R wrappers in `+st/+run/` follow the same pattern:
1. Call `commoncheck_R(rscriptdir)` — validates R path (stored in scGEAToolbox pref `scgeatoolbox/rexecutablepath`), changes to the `external/R_<toolname>/` directory, verifies packages via `require.R`
2. Write input files (`input.mat`, `positions.csv`) to the working directory
3. Run `script.R` via `pkg.RunRcode()`
4. Read output files (`.h5`, `.csv`) back into MATLAB

## Release Preparation (`copyDepsWithStructure.m`)

Run from `D:/claude_dev/` to prepare `scGEAToolbox_prj/toolbox/` for MATLAB toolbox packaging. See the detailed breakdown in prior conversation context. After running, open `scGEAToolbox_prj/scGEAToolbox.prj` in MATLAB to build the `.mltbx`.

> **Note:** This script is for scGEAToolbox release, not stGEAToolbox — stGEAToolbox is distributed directly from its source folder.

## MATLAB Version

Requires: MATLAB 2025b (same as scGEAToolbox), Statistics and Machine Learning Toolbox.
