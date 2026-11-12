# Data Availability — scatac-calib v1.0.0

This document lists every input behind the wheat benchmark, where it comes from,
which reference/coordinate system it uses, and which derived files are needed to
reproduce Table 2/3.  No raw sequencing data or full matrices are redistributed
with the package; the package ships only the 100 × 80 real-data fixture.

## 1. Wheat discovery data (primary)

| Item | Value |
|---|---|
| Biological material | *Triticum aestivum* cv. AK58, root |
| Original study | Zhang et al., *Genome Biology* 2023, **24:65**, doi:10.1186/s13059-023-02908-x (PMC10074895) |
| Repository | Genome Sequence Archive (GSA), NGDC: **CRA008788** |
| Experiment accession | **CRX537336** (“AK58-snATAC”) |
| Project accession | PRJCA012935 |
| Assay | 10x Chromium scATAC-seq (Chip E), ~10,000 nuclei; 23,518 barcodes deposited |
| Upstream in original study | Cell Ranger ATAC v2.0 + ArchR v1.0.1; 22 cell types |
| Access condition | GSA/NGDC; follow the repository terms. Raw reads require a GSA account/browser session (scripted access can hit a SecurityCheck page) |

### Processing performed here (wheat)

| Step | Result / file |
|---|---|
| QC | cells with ≥1,000 in-peak counts and ≥500 detected peaks → **8,661 cells** |
| Full peak catalogue | MACS2 v2.2.9.1, **403,519 peaks** (`results/full_alignment/AK58_peaks_peaks.narrowPeak`) |
| Alignment | bowtie2 v2.5.5, 50-bp reads, **no MAPQ filter** (deliberate; see Methods §6.1), duplicates removed with samtools v1.21 `rmdup` |
| BAM | `results/full_alignment/AK58_sorted_dedup.bam` |
| Marker peak set | **31,741 peaks** (`dataset_v3`; A 9,976 / B 10,457 / D 11,308) selected with Leiden-cluster labels (marker union) |
| Matrix | peak × cell, binary; train/val/test = 23,268 / 3,808 / 4,665 peaks × 8,661 cells; split by chromosome (chr6 = val, chr3 = test) |
| Sequence | 1,344 bp centred on each peak (±672 bp), int8 (0=A, 1=C, 2=G, 3=T), N→A |
| Reference | IWGSC RefSeq **v1.0** (the directory named “v2.1” is v1.0 by content; verified against the v1.0 `.fai`) |
| Annotation | IWGSC v1.0 HighConf GFF; chromosome names normalised (`chr1A` → `1A`) |
| Cluster labels | Leiden, 14 clusters (`dataset_v3/cell_groups.txt`) |

### Derived files needed for Table 2/3

| File | Location (server) | Content |
|---|---|---|
| `m_test.npz` | `02_results/wheat/dataset_v3/m_test.npz` | 4,665 peaks × 8,661 cells, CSR |
| `m_train.npz` | `02_results/wheat/dataset_v3/m_train.npz` | 23,268 peaks × 8,661 cells, CSR |
| `m_val.npz` | `02_results/wheat/dataset_v3/m_val.npz` | 3,808 peaks × 8,661 cells, CSR |
| `cell_groups.txt` | `02_results/wheat/dataset_v3/cell_groups.txt` | `barcode<TAB>cluster`, 8,661 rows |
| `peaks.bed` | `02_results/wheat/dataset_v3/peaks.bed` | 31,741 rows, v1.0 coordinates |
| `dataset_report.json` | `02_results/wheat/dataset_v3/dataset_report.json` | peak rule, split rule, positive rate, subgenome counts |
| `m_test.npz` etc. | `02_results/wheat/dataset_v3/` | **not redistributed**; regenerate with the scripts below |
| `ad.h5ad` | `02_results/wheat/dataset_v3/ad.h5ad` | AnnData view for the 31,741 peaks (binarised labels; sequences in `*_seqs.h5`) |

Dataset construction script: `共享/代码/脚本/build_dataset_v3.py` (in the local
project tree; not shipped in the Python package).  It documents the marker rule
(`ArchR addMarkerFeatures` style: Poisson test, log2FC ≥ 1.0, p ≤ 1e-6,
detection ≥ 50) and the split rule.  Known limitation: the original peak-calling
command/parameters were not preserved, and the marker selection used cluster
labels computed on the same 8,661 cells (the test set is enriched, not
genome-wide).

### Model predictions (needed to reproduce the numbers)

| File | Location | Content |
|---|---|---|
| Peak-selected predictions | `02_results/wheat/p1_benchmark/bootstrap/seed{0,1,2}/preds.npz` | 4,665 × 8,661 float32, key `preds` |
| Cell-selected predictions | `02_results/wheat/p1_benchmark/ci_completion/seed{0,1,2}/preds_cell.npz` | 4,665 × 8,661 float32 |
| Checkpoint metrics | `02_results/wheat/p1_benchmark/ablation/checkpoint_selection.tsv` | 3 seeds × 2 rules |
| Trained weights | `02_results/wheat/p1_benchmark/seeds/seed{0,1,2}/best_model*.pth` | scBasset checkpoints; not redistributed |

### Results referenced from this work order

| File | Location |
|---|---|
| `calib_report.json` (v1.0 audit, 3 seeds) | `02_results/wheat/p1_benchmark/calib_formal/seed{0,1,2}/` |
| `B_summary.json` / `B_summary.md` | `02_results/wheat/p1_benchmark/` |
| `bootstrap_summary.json` | `02_results/wheat/p1_benchmark/` |
| Random-peak control (3 seeds) | `02_results/wheat/p1_benchmark/random_peaks_true/seed{0,1,2}/` |
| Same-peak-set random null | `02_results/wheat/p1_benchmark/random_peaks_true/zero_model.json` |
| Shuffle control | `02_results/wheat/p1_benchmark/ablation/shuffle_control.json` |
| MAPQ / subgenome sensitivity | `02_results/wheat/dataset_v3/subgenome_effect.json`, `subgenome_effect_nomapq.json`, `dataset_v3/mapq/q0.tsv`, `q30.tsv` |
| M1 split-half reliability | `项目1/结果/minloop/m1_reliability.json` (script `项目1/代码/minloop_m1.py`) |
| 100 × 80 fixture (shipped) | `examples/real_fixture/` in this package |

## 2. Rice comparator

### 2.1 GSE214131 (main comparator, read-level only)

| Item | Value |
|---|---|
| GEO / SRA | **GSE214131** / **PRJNA883990** |
| Samples | GSM6599888 (Nip_Normal_Rep1), GSM6599889 (Nip_Normal_Rep2) |
| Layout | 2 biological replicates × 4 lanes = **8 runs**, 2 × 50 bp, ~57 GB |
| Processed h5 | `*_raw_peak_bc_matrix.h5` (177 MB each; 66,559 peaks), same samples, Cell Ranger upstream |
| Reference | **MSU7** (`osa1_r7.asm.chrs.fa.gz`; chr1 max 43,475,997 bp, which exceeds IRGSP-1.0, confirming MSU7) |
| Local paths | `项目1/数据/水稻/GSE214131/raw_fastq/GSM.../`, `项目1/数据/水稻/GSE214131/normal_rep{1,2}/` |
| Critical limitation (F-014) | Deposited FASTQ **do not contain a separate barcode/UMI read** (2 × 50 bp only). Raw reads therefore cannot rebuild a cell × peak matrix and are used **only for read-level measurements** (MAPQ, unique-alignment fraction, in-peak fragment fraction). Rice single-cell analyses use the processed h5 as an explicitly non-same-pipeline track. |
| MAPQ / read-level | `02_results/rice/read_level/` (`*.mapq.tsv`, incl. `WHEAT_AK58.mapq.tsv` for the wheat sample) |

### 2.2 GSE252040 (backup; not used for the main comparison)

| Item | Value |
|---|---|
| GEO / SRA | **GSE252040** / **PRJNA1007577** |
| Layout | 4 runs, 2 × 150 bp, ~256 GB |
| Reference | **IRGSP-1.0** (Ensembl Plants release-57) |
| Status | paused; backup only. Same missing-barcode limitation (F-014) |

> **Never mix MSU7 and IRGSP-1.0 coordinates in one comparison.**  The main
> rice track uses MSU7; IRGSP-1.0 is only for the unused backup series.

## 3. Annotation and reliability inputs

| Data | Content | Source |
|---|---|---|
| Zhang Table S2 | 22 clusters × top-100 marker genes | Zhang et al. 2023 supplementary |
| Zhang Table S8 | cluster × A/B/D gene × RNA expression | Zhang et al. 2023 supplementary; used for the ATAC–RNA exploratory check |
| Homoeolog triads | 18,407 A/B/D triads (`triad` table + 109 MB windows) | derived from the v1.0 annotation |
| IWGSC ID correspondence | `iwgsc_refseq_all_correspondances.csv` (26 MB; v1.0/v1.1/v2.1) | IWGSC |
| Wheat snRNA (optional) | GSA CRA008788 / CRX537334 | Zhang et al. 2023; not used for the main claims |

The v1.1 → v1.0 gene-ID mapping is required because Table S2/S8 use v1.1 IDs
(`TraesCS...02G...`) while our peaks/GFF use v1.0 (`...01G...`); direct
intersection is zero.  1,736/1,757 Table S2 markers (98.8%) map successfully.

## 4. What is and is not redistributed

- **Redistributed in the Python package:** the 100-peak × 80-cell real-data
  fixture (`examples/real_fixture/`: labels, 80 cell columns of the train
  matrix, seed-0 predictions, groups, BED, selection manifest) and the fully
  synthetic demo.
- **Not redistributed:** raw FASTQ/BAM, the full 4,665 × 8,661 and
  23,268 × 8,661 matrices, sequences, model weights, Zhang Table S2/S8,
  rice raw reads, and the rice processed h5.
- **Access:** wheat raw data via GSA CRA008788 under the original repository
  terms; rice raw/processed data via GEO GSE214131 / SRA PRJNA883990 under the
  original repository terms.  Redistribution of derived matrices must respect
  those terms.

## 5. Suggested Data Availability statement (fill placeholders)

> The wheat scATAC-seq data analysed here are from Zhang et al. (2023,
> *Genome Biology* 24:65; GSA CRA008788 / CRX537336, project PRJCA012935) and
> were processed as described in Methods.  The rice comparator is GEO
> GSE214131 / SRA PRJNA883990; its deposited FASTQ lack a separate barcode read,
> so raw reads are used only for read-level statistics and single-cell analyses
> use the deposited processed h5 (MSU7 coordinates).  Derived matrices, model
> predictions and per-seed audit outputs are available at `{DATA_REPOSITORY_URL}`
> (`{DATA_DOI}`), and the 100 × 80 real-data fixture plus expected outputs are
> included in the `scatac-calib` v1.0.0 archive (`{ZENODO_DOI}`).

Do not paste until `{DATA_REPOSITORY_URL}` and `{DATA_DOI}` are decided.
