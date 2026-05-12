# bioepy — OTU Clustering & Genomic Analysis Library
**Author:** Epicoco Andrea  
**Status:** Development / Testing

---

## Overview

`bioepy` is a Python library for bioinformatics analysis, currently focused on two independent workflows:

1. **OTU Clustering** — Reference-based clustering of Paired-End sequencing reads using a k-mer Jaccard similarity approach
2. **Gene/Transcript Analysis** — Extraction of gene sequences from genomic assemblies with transcription and translation support

The library is contained in `bioepy.py` and used via two entry-point scripts: `main.py` (OTU clustering pipeline) and `class_test.py` (gene/transcript testing).

---

## Project Structure

```
otu_cluster/
├── main.py                        # Entry point: OTU clustering pipeline
├── class_test.py                  # Entry point: Gene/Transcript analysis test
├── bioepy.py                      # Core library (all classes and functions)
├── file/
│   ├── merging_test1.fq.gz        # PE reads - forward (R1) — test input
│   ├── merging_test2.fq.gz        # PE reads - reverse (R2) — test input
│   ├── ref.fa                     # Reference sequences for OTU clustering
│   ├── chr7.fa                    # Chromosome 7 sequence (T2T assembly)
│   ├── ensembl.txt.gz             # Ensembl annotation table (for gene extraction)
│   └── caratteristiche_sequenze_merging.txt  # Notes on test files
├── cluster_file_k8/               # Output: clustering results for k-mer length 8
│   └── OTU_clustered_seq_*.txt
└── old_version/                   # Archived previous implementations
```

---

## Requirements

- Python 3.10+
- `pandas`
- `dask`
- `regex`

Install dependencies:
```bash
pip install pandas dask regex
```

---

## Workflow 1 — OTU Clustering (`main.py`)

### What it does

The pipeline performs **reference-based OTU (Operational Taxonomic Unit) clustering** on merged Paired-End reads in three steps:

```
PE reads (R1 + R2)
        │
        ▼
1. Quality filtering & FastQ parsing
        │
        ▼
2. Merging of Paired-End reads
        │
        ▼
3. OTU clustering vs. reference sequences (Jaccard similarity)
        │
        ▼
Output: cluster_file_k{N}/OTU_clustered_seq_*.txt
```

### Configuration

Before running, edit `main.py` and set `FOLD_PATH` to your project folder:

```python
FOLD_PATH = "C:/path/to/your/project/folder/"
```

The folder must contain a `file/` subfolder with:

| File | Description |
|------|-------------|
| `merging_test1.fq.gz` | R1 Paired-End reads (FastQ, gzipped) |
| `merging_test2.fq.gz` | R2 Paired-End reads (FastQ, gzipped) |
| `ref.fa` | Reference sequences in FastA format |

You can also adjust the k-mer range used for clustering:

```python
KMER_LEN_MIN = 8   # minimum k-mer length
KMER_LEN_MAX = 12  # maximum k-mer length
```

### Run

```bash
python main.py
```

---

### Step 1 — FastQ Parsing (`ParsingFastq`)

Each FastQ file is parsed in parallel using multiprocessing. For each read the parser:

- Validates the sequence ID (must start with `@`)
- Validates the DNA sequence (only `ATCGN` characters)
- Validates the quality score length
- Auto-detects the Phred encoding: **Phred-33** (Sanger/Illumina ≥1.8) or **Phred-64** (Illumina ≤1.7)

After parsing, sequences are filtered by **Expected Error (EE)**:

- Default threshold: `EE ≤ 2.0`
- If fewer than 85% of sequences pass, the user is prompted to increase the threshold (`+0.5` per step)
- Duplicate sequences are removed, keeping the one with the highest average quality score

**Instance attributes available after parsing:**

| Attribute | Description |
|-----------|-------------|
| `ids` | List of sequence identifiers |
| `seqs` | List of DNA sequences |
| `qss` | List of quality score strings |
| `qss_num` | Quality scores converted to numeric (Phred) |
| `exp_errs` | Binary list: `1` = valid, `0` = discarded |
| `valid_seqs` | Filtered unique sequences `[[id, seq], ...]` |

---

### Step 2 — PE Merging (`MergingPE`)

Paired-End reads are merged by finding the **overlap** between the forward read (R1) and the reverse complement of the reverse read (R2).

- Only reads sharing the same accession ID are paired
- Overlap is searched iteratively from the full length down to a minimum of `min_overlap=8` bases
- Up to `max_err=3` mismatches are tolerated in the overlap region (fuzzy matching via `regex`)
- The merged sequence is: `seq1[:-overlap] + seq2[overlap_start:]`

---

### Step 3 — OTU Clustering (`ClusterOTU`)

Each merged sequence is compared to every reference sequence using the **Jaccard similarity index** based on k-mers.

Two methods are used depending on the reference sequence:

| Method | Condition | Description |
|--------|-----------|-------------|
| **Standard Jaccard** (`S`) | No `N` bases in reference | `|intersection(k1, k2)| / |k1|` |
| **Weighted Jaccard** (`W`) | `N` bases present in reference | k-mers weighted by `alpha^n_count` (penalizes ambiguous bases) |

A merged sequence is assigned to a reference cluster if the similarity index is **≥ 0.97**.

The clustering runs over all k-mer lengths from `KMER_LEN_MAX` down to `KMER_LEN_MIN`, and is **fully parallelized** using Python `multiprocessing.Pool` (uses all CPUs − 1).

### Output

For each k-mer length, a folder `cluster_file_k{N}/` is created containing one file per CPU worker:

```
cluster_file_k8/
├── OTU_clustered_seq_0.txt
├── OTU_clustered_seq_1.txt
└── ...
```

Each output file lists, for each reference sequence that attracted at least one merged read:

```
>reference_sequence_id
[Wji=0.975] [merged_read_id]    TGGGGA...AACA
[Sji=0.981] [merged_read_id]    ATCGGA...TTCA
```

- `W` / `S` = method used (Weighted / Standard Jaccard)
- `ji` = Jaccard index value
- Files with no matches (< 5 bytes) are automatically deleted

---

## Workflow 2 — Gene & Transcript Analysis (`class_test.py`)

> ⚠️ **Currently in testing.** The script is hardcoded to extract and analyze the **CFTR gene** from chromosome 7.

### What it does

```
chr7.fa (chromosome sequence)  +  ensembl.txt.gz (annotation)
        │
        ▼
Gene: extract annotated regions → sequences from the chromosome
        │
        ▼
Transcript: DNA → mRNA (3 reading frames)
        │
        ▼
Translation: mRNA → protein (only for protein_coding genes)
```

### Classes

#### `ParsingFasta`
Parses a FastA file using Dask. Populates:
- `ids` — list of sequence headers (without `>`)
- `seqs` — list of corresponding sequences

#### `Gene`
Extracts all annotated features for a given gene from an Ensembl-format annotation table, matched against the chromosome sequence.

- Assembly: **T2T CHM13.v2.0** (configurable via class methods)
- Produces `gene_info_df`: a DataFrame with annotation metadata and the extracted DNA sequence for each feature

#### `Transcript(Gene)`
Extends `Gene` with:
- `transcription()`: generates mRNA in all **3 reading frames** (`transcript_seq0/1/2`)
- `translation()`: translates to protein for `protein_coding` genes (`protein_seq0/1/2`)

### Run

```bash
python class_test.py
```

Requires `chr7.fa` and `ensembl.txt.gz` in the configured `FILE_FOLD_PATH`.

---

## Current Limitations & Known Issues

- `FOLD_PATH` and `FILE_FOLD_PATH` are **hardcoded** in `main.py` and `class_test.py` — must be edited before use. The interactive file selection mode exists in the code but is currently commented out.
- `cpu.parse.py` is present but **empty** — placeholder for a future CPU parsing module.
- The `class_test.py` is fixed to the CFTR gene (`ENSG00000001626`) and chromosome 7 — not yet generalized.
- `old_version/` contains previous iterations of the clustering algorithm — kept for reference, not maintained.
- The weighted Jaccard index has a commented-out alternative implementation in `bioepy.py` — the current version normalizes over `k1` only (asymmetric similarity).
