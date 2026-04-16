# MixOmics — Transcriptomic & Metabolomic Statistical Analysis
**Author:** Epicoco Andrea  
**Version:** 1.0

---

## Overview

**MixOmics** is a command-line tool for the integrated statistical analysis of **transcriptomic** and **metabolomic** data generated from RESOLUTE experimental results.

The tool processes files where each gene or metabolite is associated with a **p-value** and **q-value** (adjusted p-value) that quantify the statistical significance of expression/abundance variation between a **control sample** and a **treated/perturbed sample** (e.g. Knock-Out vs Over-Expression).

---

## Requirements

- **Windows** operating system
- **Conda / Miniconda** (automatically installed if not found)
- Internet connection (only required on first run, for Miniconda or environment setup)

> On first launch, the script will automatically:
> 1. Check for an existing Conda installation
> 2. Download and install Miniconda if not found
> 3. Create the `mixomics` conda environment from `environment.yml`

---

## Input Data

The tool expects two CSV files per experiment, placed in a subfolder of `./data/` named according to the following convention:

```
./data/{gene_name}_{cell_line}_{exp_type}/
```

| File | Description |
|------|-------------|
| `{gene_name}_gene_all.csv` | Transcriptomic data (DESeq2 output with p-value and q-value per gene) |
| `{gene_name}_mtbl_all.csv` | Metabolomic data (p-value and q-value per metabolite) |

Each file must contain, for every gene or metabolite, the **p-value** and **q-value** expressing the significance of variation between the control and the treated condition.

### Example folder structure

```
./data/
└── slc25a36_hct116_kooe/
    ├── slc25a36_gene_all.csv
    └── slc25a36_mtbl_all.csv
```

> **Note:** If the folder does not exist yet, the program will create it automatically and prompt you to drop the input files before proceeding.

---

## Usage

```
mixomics.bat -gn <gene_name> [-cl <cell_line>] [-et <exp_type>]
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `-gn`, `--gene_name` | Name of the knocked-out (KO) gene |

### Optional arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `-cl`, `--cell_line` | Name of the cell line used | `hct116` |
| `-et`, `--exp_type` | Type of experiment comparison | `kooe` |
| `-h`, `--help` | Show help message | — |

### Experiment types (`-et`)

| Code | Comparison |
|------|------------|
| `kooe` | KO vs Over-Expression (OE) |
| `koct` | KO vs Control |
| `wtko` | Wild-Type (WT) vs KO |

---

## Examples

```bat
:: Minimal usage (defaults: hct116, kooe)
mixomics.bat -gn slc25a36

:: Full explicit call
mixomics.bat -gn slc25a36 -cl hct116 -et kooe

:: Using long argument names
mixomics.bat --gene_name slc25a36 --cell_line hct116 --exp_type kooe
```

---

## How It Works

```
Launch mixomics.bat
        │
        ▼
Check / install Conda & environment
        │
        ▼
Read input CSV files from ./data/{gene_name}_{cell_line}_{exp_type}/
        │
        ├── {gene_name}_gene_all.csv   ← transcriptomic (DESeq2: p-value, q-value per gene)
        └── {gene_name}_mtbl_all.csv   ← metabolomic (p-value, q-value per metabolite)
        │
        ▼
Statistical analysis
  - Significance filtering by p-value / q-value
  - Integration of transcriptomic and metabolomic layers
        │
        ▼
Output results
```

---

## Project Structure

```
mixomics/
├── mixomics.bat          # Launcher script (entry point)
├── environment.yml       # Conda environment definition
├── dist/
│   └── main.py           # Main Python analysis script
└── data/
    └── {gene}_{line}_{type}/
        ├── {gene}_gene_all.csv
        └── {gene}_mtbl_all.csv
```

---

## Output Files & Interpretation

All output files are saved in the same experiment folder: `./data/{gene}_{cell_line}_{exp_type}/`

---

### Transcriptomic Results

| File | Description |
|------|-------------|
| `{gene}_gene_de.csv` | All differentially expressed (DE) genes |
| `{gene}_gene_down.csv` | DE genes that are **downregulated** in the treated condition |
| `{gene}_gene_up.csv` | DE genes that are **upregulated** in the treated condition |
| `{gene}_gene_marker.csv` | DE genes with **\|log2FC\| > 2** — strong expression changes, candidate markers |
| `{gene}_kegg_analysis.log` | Full execution log — useful for tracing analysis steps and getting a quick overview of results without opening individual files |

---

### Metabolomic Results

| File | Description |
|------|-------------|
| `{gene}_mtbl_da.csv` | All differentially abundant (DA) metabolites |
| `{gene}_mtbl_down.csv` | DA metabolites that are **downregulated** |
| `{gene}_mtbl_up.csv` | DA metabolites that are **upregulated** |
| `{gene}_mtbl_marker.csv` | DA metabolites with **\|log2FC\| > 2** *(not always generated)* |

---

### KEGG Pathway Analysis

#### Summary tables

These TSV files report KEGG pathways enriched by DE genes (up/down) or DA metabolites, and contain the following columns:

| Column | Description |
|--------|-------------|
| `kpid` | KEGG Pathway ID |
| `kpname` | KEGG Pathway name |
| `nmp` | Number of metabolites in the pathway |
| `nmda` | Number of DA metabolites found in the pathway |
| `ngp` | Number of genes in the pathway |
| `ngde` | Number of DE genes found in the pathway |
| `odd` | Odds ratio from the contingency table |
| `pvalue` | Statistical significance of the difference between control and treated |
| `qvalue` | Adjusted p-value (corrected for multiple testing) |
| `kgid` | List of KEGG gene IDs in the pathway |
| `kgsymbol` | List of KEGG gene symbols |
| `kgde` | KEGG IDs of DE genes present in the pathway |
| `kmid` | KEGG compound IDs of metabolites in the pathway |
| `kmname` | Names of metabolites in the pathway |
| `kmda` | KEGG IDs of DA metabolites in the pathway |

| File | Content |
|------|---------|
| `{gene}_kegg_path_info_gene_down.tsv` | KEGG pathways statistically significant for **downregulated** genes (Gene Ontology-driven) |
| `{gene}_kegg_path_info_gene_up.tsv` | KEGG pathways statistically significant for **upregulated** genes (Gene Ontology-driven) |
| `{gene}_kegg_path_info_mtbl.tsv` | All KEGG pathways containing **at least one DA metabolite** from the metabolomic analysis |

---

### Gene Ontology Results

| Folder | Description |
|--------|-------------|
| `{gene}_go_results/` | Results of Gene Ontology enrichment analysis, performed via **g:Profiler** |

---

### KEGG Enrichment Folders (for pathway visualization)

These folders contain query files to visualize DE genes and DA metabolites directly on KEGG pathway graphs.

| Folder | Content |
|--------|---------|
| `{gene}_kegg_enrich_down/` | Query files for pathways identified as **downregulated** by GO analysis |
| `{gene}_kegg_enrich_up/` | Query files for pathways identified as **upregulated** by GO analysis |
| `{gene}_kegg_enrich_mtbl/` | Query files for pathways identified through **DA metabolites** |

Each folder contains pairs of files per pathway:

| File | Description |
|------|-------------|
| `{pathway_name}_{kegg_id}.ge` | Genes listed by their KEGG gene ID, each associated with its log2FC value |
| `{pathway_name}_{kegg_id}.me` | Metabolites listed by their KEGG compound ID, each associated with its log2FC value |

#### How to visualize on KEGG pathway graphs

1. Open the KEGG pathway page for the pathway of interest
2. In the left side panel, click the **"+"** icon next to the **"color"** label
3. Load the `.ge` file — tick the option **"Numerical values converted to color gradation"** and choose your preferred color scale
4. Repeat step 3 for the `.me` file
5. Genes and metabolites will be colored according to their log2FC, allowing immediate visual identification of up- and downregulated elements within the pathway

---

### Over-Representation Analysis (ORA)

| Folder | Description |
|--------|-------------|
| `{gene}_ora/` | Results of Over-Representation Analysis on KEGG pathways. Identifies which pathways are significantly enriched by DA metabolites, based purely on **quantitative overlap** (not on the physiological effect of individual metabolite changes). Useful as a complementary, statistics-driven perspective to the pathway enrichment analysis. |

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `FileNotFoundError` on CSV | Check that the input files exist in the correct folder and are named exactly as expected |
| `conda not found` | Miniconda will be downloaded automatically; ensure internet access on first run |
| Environment creation fails | Verify that `environment.yml` is present in the same folder as `mixomics.bat` |
| Input with trailing spaces | The launcher automatically trims spaces from all input arguments |
