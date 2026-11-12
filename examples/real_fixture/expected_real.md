# Real-data fixture — expected output (v1.0.0)

This directory is a **real-data fixture**, not the synthetic demo.  It is a
fixed-index subset of the wheat `dataset_v3` test split, taken read-only from
`L-wdl:/root/wheat_project/02_results/wheat/dataset_v3` and the seed-0
peak-selected model predictions
(`p1_benchmark/bootstrap/seed0/preds.npz`).  It exists so that a third party
can verify the CLI end to end on real scATAC-seq data without downloading the
full 4,665 × 8,661 matrices.

> This fixture is an interface/reproduction test.  Its AUROCs are computed on
> 100 peaks × 80 cells and **must not be quoted as dataset-level performance**.
> The paper's numbers come from the full test split (`docs/reproducibility.md`).

## Inputs

| File | Shape / content | SHA256 |
|---|---|---|
| `m_test_fixture.npz` | test labels, CSR, **100 peaks × 80 cells**, binary, 2.1875% accessible | `55c23785…14a10` |
| `m_train_fixture.npz` | training labels, CSR, **23,268 peaks × 80 cells** (all train peaks, selected cells only) | `99051b9d…89995` |
| `preds_fixture.npz` | key `preds`, float32, **100 × 80**, seed-0 peak-selected model | `4632c251…67c24` |
| `groups_fixture.txt` | 80 rows `barcode<TAB>cluster`, 14 Leiden clusters | `dfce0db6…996b9` |
| `peaks_fixture.bed` | 100 rows `chr,start,end,name` (chr3; 1,344-bp window centred on the peak) | `43cd66bf…40712` |

Full hashes, exact row indices and the selection rule are in
`selection_manifest.json`.

### How the subset was fixed

- **Peaks:** `np.unique(np.round(np.linspace(0, 4664, 100)))` over the 4,665
  test peaks (chr3), i.e. every ~47th test peak.  Row order in `m_test_fixture`
  is this ascending index order; `peaks_fixture.bed` is aligned row by row.
- **Cells:** the first 5 cells of each of the 14 Leiden clusters (the first 10
  clusters in numeric order get 6), then sorted by original cell index.
  All 14 groups are represented: 10 × 6 + 4 × 5 = 80.
- **Predictions:** the same peak and cell indices from seed-0 peak-selected
  `preds.npz`.
- **Training matrix:** all 23,268 training peaks, restricted to the same 80
  cells.  Train and test peaks are disjoint by chromosome, so
  `zero_model.cell_freq_source = "train"` and the oracle is
  `test(label-informed)` (an upper bound, not a trainable result).

## Command

```bash
cd examples/real_fixture
./run_real_fixture.sh            # writes fixture_out/
```

Equivalent explicit form:

```bash
scatac-calib audit \
  --matrix m_test_fixture.npz \
  --preds preds_fixture.npz \
  --train-matrix m_train_fixture.npz \
  --groups groups_fixture.txt \
  --out fixture_out \
  --n-boot 1000 --depth-bins 5 --seed 20261112
```

## Reference numbers

Generated with **scatac-calib 1.0.0**, Python 3.14.7, numpy 2.5.4,
scipy 1.18.1, scikit-learn 1.9.1.

| Quantity | Per-peak | Per-cell |
|---|---|---|
| Model AUROC (95% CI, mean bootstrap) | **0.747824** [0.703130, 0.789679] | **0.694087** [0.645265, 0.742089] |
| Frequency null (zero model) | 0.683415 | 0.808941 |
| ΔAUROC (model − null) | **+0.064409** [+0.028490, +0.100652], p < 0.001 | **−0.114853** [−0.170263, −0.060559], p < 0.001 |
| Group-level oracle (source `test(label-informed)`) | 0.934845 (14 groups) | 0.955008 (14 groups) |
| Depth-matched per-peak AUROC (5 bins) | 0.623545 | — |

Other fields:

- `ap` = **0.077476**
- `metrics.n_valid_peaks` = 73 of 100 (27 peaks are all-zero or all-one in this
  subset and are skipped)
- `metrics.n_valid_cells` = 64 of 80 (16 cells have no positive peak among the
  100 and are skipped)
- `zero_model.cell_freq_source` = `"train"`
- `oracle.missing_cells` = 0
- `depth_matched.n_bins` = 5 (4 bins contributed; one rank-bin contained no
  peak with both classes)
- `depth_matched.per_bin_auc` =
  [0.583471, 0.593786, 0.628810, 0.688112]
- `depth_matched.cell_depth_median` = 2.0

## Tolerance

Point estimates (model / zero / oracle / depth AUROC, AP) should match to
**1e-6** on a given numpy release; bootstrap CIs and p-values to **1e-3**.
`p < 0.001` prints as `0.0` because p is clipped at 0 when no bootstrap
replicate crosses zero.

## Reading the numbers (why the per-cell Δ is negative here)

Do **not** interpret this subset as a failure or a result.  Three properties of
the 100 × 80 slice explain it:

1. The per-cell null scores peaks by their detection frequency **in the same
   100-peak test matrix**.  With only 100 peaks, each label contributes a large
   share (1/100) to the null score; this transductive self-scoring is strongest
   at small peak counts.  On the full 4,665-peak test set the same null is
   0.6878 and the model is 0.7048 (Δ +0.0168).
2. The 80 cells are very sparse (median depth 2 across the 100 selected peaks),
   16 of them have no positive label and are skipped.
3. The predictions are a subset of a model trained and selected on the full
   4,665-peak set; ranking only every ~47th peak is not the task the model was
   selected for.

The fixture's job is to prove that installation, matrix parsing, null, oracle,
delta/bootstrap and depth stratification all run on genuine data and reproduce
to the stated tolerance.  Scientific interpretation belongs to the full-split
run.

## Difference from the synthetic demo

| | `examples/demo/` | `examples/real_fixture/` |
|---|---|---|
| Data | synthetic (seeded generator) | real wheat `dataset_v3` subset |
| Size | 100 peaks × 80 cells | 100 peaks × 80 cells |
| Train matrix | synthetic, same generative process | real, 23,268 disjoint train peaks |
| Groups | 4 synthetic groups × 20 cells | 14 real Leiden clusters |
| Purpose | self-contained tutorial, no data access needed | provenance/reproduction check on real data |
| Output | `demo_out/` | `fixture_out/` |
