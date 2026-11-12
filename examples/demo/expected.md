# Demo expected output (v1.0.0)

A third party can use this file to check that their installation reproduces
the reference run.  Numbers were generated with **scatac-calib 1.0.0**, the
default `--seed 20261112`, `--n-boot 200`, `--depth-bins 5`.

## Inputs

| File | Shape / content |
|---|---|
| `m_demo.npz` | test matrix, CSR, **100 peaks × 80 cells**, binary, 14.61% accessible |
| `m_train_demo.npz` | independent training matrix, same shape and generative process |
| `preds_demo.npz` | key `preds`, float32, true rate + N(0, 0.10) noise |
| `groups_demo.txt` | 80 rows `cell_<id><TAB>G1..G4`, 20 cells per group |

The data are synthetic: 4 groups, 25 group-marker peaks each (13/13/12/12 for
G1–G4).  Regenerate deterministically with:

```bash
python examples/demo/make_demo_data.py --seed 20261112
```

## Command

```bash
cd examples/demo
./run_demo.sh            # writes demo_out/
```

Equivalent explicit form:

```bash
scatac-calib audit \
  --matrix m_demo.npz \
  --preds preds_demo.npz \
  --train-matrix m_train_demo.npz \
  --groups groups_demo.txt \
  --out demo_out \
  --n-boot 200 --depth-bins 5
```

## Reference numbers

| Quantity | Per-peak | Per-cell |
|---|---|---|
| Model AUROC (95% CI, mean bootstrap) | **0.62766** [0.58561, 0.66816] | **0.74127** [0.72379, 0.75980] |
| Frequency null (zero model) | 0.50583 | 0.68782 |
| ΔAUROC (model − null) | **+0.12183** [+0.08141, +0.16493], p < 0.001 | **+0.05345** [+0.03605, +0.07110], p < 0.001 |
| Group-level oracle (source `train`) | 0.64660 (4 groups) | 0.75501 (4 groups) |
| Depth-matched per-peak AUROC (5 bins) | 0.66879 | — |

Other fields:

- `ap` = **0.44714**
- `zero_model.cell_freq_source` = `"train"`
- `n_valid_peaks` = 100, `n_valid_cells` = 80
- `depth_matched.per_bin_auc` =
  [0.748995, 0.663446, 0.647717, 0.630696, 0.653080]
- `depth_matched.cell_depth_median` = 14.0
- `oracle.missing_cells` = 0

## Tolerance

Bootstrap CIs and p-values vary in the last digits across numpy versions
because `np.random.default_rng` streams for `integers` are stable but BLAS
averaging may differ.  Point estimates (AUROC, zero model, oracle, depth)
should match to **1e-6** on a given numpy release; CIs to **1e-3**.
`p < 0.001` is expected to print as `0.0` because p is clipped at 0 when no
bootstrap replicate crosses zero.

## Why these numbers are meaningful

- The frequency null is ≈ 0.5 per peak (synthetic data have no depth
  gradient), so the model's +0.12 per-peak Δ is real signal in the generating
  process, not frequency.
- Per-cell null is 0.688 because cells differ in overall accessibility across
  the 25 marker peaks; this is exactly the "free" component the tool exists
  to expose.
- The group-level oracle (0.755 per cell) is above the model (0.741): knowing
  the cell group is still slightly better than the noisy predictions, which
  demonstrates the intended oracle interpretation.
