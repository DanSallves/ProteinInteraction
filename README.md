# Interaction Scores using AlphaFold3 protein prediction

This repository contains two R scripts for batch-processing **AlphaFold3** prediction output archives (`.zip`) in a local folder. The scripts are designed to help with two common post-processing tasks (1) to **calculate an interaction score** from the inter-chain **Predicted Aligned Error (PAE)** matrix and (2) to **generate PAE heatmaps** for each prediction and export them as multi-panel PNG pages. These scripts are useful for screening large sets of AlphaFold3 predictions for potential protein-protein interactions and for visually inspecting model confidence across and between chains.

---

## Overview

Both scripts work by scanning the current working directory for AlphaFold3 result archives (`*.zip`). For each archive, they:

- extract the zip into a temporary directory,
- locate the required JSON files,
- read the **PAE matrix** from the `*_full_data_0.json` file,
- read the chain sequences from the `*_job_request.json` file,
- and then either:
  - compute an interaction score, or
  - generate a PAE heatmap.

The scripts assume each prediction contains **two protein chains**, and they treat:

- **sequence 1** as the **prey**, and
- **sequence 2** as the **bait**.

---

## Scripts

### 1. Interaction score calculation

This script batch-processes zipped AlphaFold3 predictions and calculates an **interaction score** from the inter-chain region of the PAE matrix.

#### What it does

For each `.zip` file in the working directory, the script:

1. extracts the archive to a temporary directory,
2. finds:
   - `*_full_data_0.json`
   - `*_job_request.json`
3. reads the two protein sequences from the job request file,
4. determines the length of each protein,
5. extracts the PAE matrix,
6. isolates the **inter-chain PAE block**,
7. converts that block into a summed score,
8. stores the result in a summary table.

#### Output

The script writes a CSV file:

```text
SPRYonly_Contacts_OAS2.csv
```

with columns:

- `Gene` - derived from the zip filename
- `Contacts` - the calculated interaction score
- `Prey` - length of protein 1
- `Bait` - length of protein 2

---

### 2. PAE plot generation

This script batch-processes the same type of AlphaFold3 zip archives and creates **PAE heatmaps** for each prediction.

#### What it does

For each `.zip` file in the working directory, the script:

1. extracts the archive to a temporary directory,
2. finds:
   - `*_full_data_0.json`
   - `*_job_request.json`
3. reads the PAE matrix,
4. reads the lengths of the two protein sequences,
5. converts the matrix into long format,
6. generates a heatmap using `ggplot2`,
7. draws horizontal and vertical boundary lines marking the end of the first sequence,
8. stores all plots and exports them in pages of **18 plots per PNG** arranged in a **6 x 3** layout.

#### Output

The script creates PNG files such as:

```text
PAE_Page_1.png
PAE_Page_2.png
...
```

Each PNG contains up to 18 PAE plots.

---

## Input requirements

Both scripts expect the working directory to contain one or more AlphaFold3 output zip files.

Each zip archive should include at least:

- a file ending with:

```text
_full_data_0.json
```

- and a file ending with:

```text
_job_request.json
```

The scripts search recursively inside each extracted archive.

---

## Dependencies

### Required R packages

For the interaction score script:

```r
library(jsonlite)
library(tools)
library(jsonlite)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(tools)
```

---

## Usage

### Step 1. Place input zip files in one folder

Put all AlphaFold3 `.zip` output files into a single directory.

### Step 2. Open the appropriate script and set the working directory

In the interaction score script, the working directory is explicitly set using:

```r
setwd("/path/to/your/folder")
```

In the PAE plotting script, no `setwd()` is included, so the script uses the current R working directory. You can either:

- run the script from the correct folder, or
- add your own `setwd()` line at the top.

### Step 3. Run the script in R or RStudio

For example:

```r
source("calculate_interaction_score.R")
source("plot_pae_panels.R")
```

Replace the filenames above with the actual names used in your repository.

---

## Outputs

### Interaction score script

Produces a summary CSV file similar to:

```text
SPRYonly_Contacts_OAS2.csv
```

Example output table:

| Gene | Contacts | Prey | Bait |
|------|----------|------|------|
| GENE1 | 123456.7 | 350 | 120 |
| GENE2 | 110983.5 | 350 | 98 |

### PAE plotting script

Produces one or more image files:

```text
PAE_Page_1.png
PAE_Page_2.png
...
```



