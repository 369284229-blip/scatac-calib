# Release checklist — scatac-calib v1.0.0

Prepared by worker 050 on 2026-11-12.  Everything that can be done without the
authors' decisions is done: the repository is initialised, the release commit is
tagged `v1.0.0`, the availability documents are drafted, the real-data fixture
and the Table 3 executable mapping are included, and the wheel/sdist are
rebuilt.  The items below need **you** (the authors), because they involve
identity, hosting and permanent identifiers.

Nothing in this file has been faked: every `{TBD}` is still a placeholder.

## 1. Decisions only the authors can make

| # | Decision | Where it is used | Suggested default |
|---|---|---|---|
| 1 | **Public or private repository?** | `CODE_AVAILABILITY.md`, manuscript Code Availability | Public (reviewers must be able to run it) |
| 2 | **Host and repository URL** | `{REPOSITORY_URL}` in `CODE_AVAILABILITY.md` | GitHub, e.g. `https://github.com/<org>/scatac-calib` |
| 3 | **Zenodo DOI** (or other archive) | `{ZENODO_DOI}` in both availability docs | Enable the GitHub–Zenodo integration, then create a release for `v1.0.0`; the DOI is minted automatically |
| 4 | **Author list and affiliations** | `CITATION.cff`, `pyproject.toml` authors, LICENSE copyright | Replace “Wheat scATAC Benchmark Team” with the real authors |
| 5 | **License holder / copyright line** | `LICENSE` | Institution or the corresponding author; MIT is already chosen |
| 6 | **Contact e-mail** | `CODE_AVAILABILITY.md` | Corresponding author |
| 7 | **Data repository URL/DOI** (for derived matrices/predictions) | `DATA_AVAILABILITY.md` `{DATA_REPOSITORY_URL}` / `{DATA_DOI}` | Zenodo/figshare, or state “available on request” if the source terms forbid redistribution |
| 8 | **Whether derived matrices can be redistributed** | `DATA_AVAILABILITY.md` §4 | Check the GSA/GEO terms before uploading |

## 2. Commands you run once the URL/DOI exist

From `代码/scatac_calib/`:

```bash
# 1. Create the remote repository (empty; do not let GitHub add a README).
#    Example for GitHub over SSH (required on the campus network):
git remote add origin git@github.com:<org>/scatac-calib.git

# 2. Push the branch and the tag.  The tag already exists locally.
git push -u origin main          # or: git push -u origin master
git push origin v1.0.0

# 3. Archive: on Zenodo, enable the repository, then create a GitHub release
#    for tag v1.0.0.  Copy the DOI into:
#      CODE_AVAILABILITY.md  -> Archived DOI (Zenodo)
#      DATA_AVAILABILITY.md  -> {ZENODO_DOI}
#      CITATION.cff          -> (optional) identifiers

# 4. Fill the remaining {TBD}s (author list, repository URL, contact,
#    data repository URL/DOI) and commit:
git add -A
git commit -m "Fill release identifiers for v1.0.0"
git tag -f v1.0.0                # only if you want the tag to include the edits
```

> On the current campus network, GitHub HTTPS is blocked but SSH over 443 works
> (`~/.ssh/config` already maps `github.com` to `ssh.github.com:443`).  Use the
> SSH remote, not HTTPS.

If the repository must stay private, use Zenodo's “restricted” record or an
institutional archive, and state the access condition in the manuscript.

## 3. Pre-release verification (already run by worker 050; re-run after edits)

```bash
# Tests (9 tests, no pytest needed)
python -m unittest discover -s tests -v

# Synthetic demo
cd examples/demo && ./run_demo.sh && cd ../..

# Real-data fixture
cd examples/real_fixture && ./run_real_fixture.sh && cd ../..

# Package build and hashes
python -m build --no-isolation
shasum -a 256 dist/*
```

Update `dist/SHA256SUMS` (or the table in `回执_050.md`) if you rebuild.

## 4. Content checks before pushing

- [ ] `COMMIT` contains the output of `git rev-parse HEAD` for the tagged commit.
- [ ] `VERSION` = `1.0.0`; `scatac_calib.__version__` = `1.0.0`;
      `pyproject.toml` version = `1.0.0`; tag `v1.0.0` exists.
- [ ] `CODE_AVAILABILITY.md` and `DATA_AVAILABILITY.md` have no remaining `{TBD}`
      **except** the ones you deliberately leave for the manuscript.
- [ ] No secrets, credentials, `.env`, SSH keys, or absolute personal paths in
      the tree (`grep -rn "/Users/\|/root/\|password\|token" .` should return
      only intentional documentation examples).
- [ ] The fixture is the intended fixed subset: `examples/real_fixture/selection_manifest.json`
      lists 100 peak indices and 80 cell indices and their SHA256.
- [ ] The synthetic demo is still present (`examples/demo/`), not replaced by the
      fixture.
- [ ] `dist/scatac_calib-1.0.0.tar.gz` and `.whl` contain the new files
      (`CODE_AVAILABILITY.md`, `DATA_AVAILABILITY.md`, `RELEASE_CHECKLIST.md`,
      `COMMIT`, `VERSION`, `docs/table3_executable.md`,
      `examples/real_fixture/*`).
- [ ] Manuscript-side fixes (not done here, to respect the “no manuscript edits”
      rule): Methods §6.10 says v0.1 and commit `{TBD}`; Fig 1D/Fig 6/README say
      “items 1–8” while Table 3 has 15 rows; Table 2 primary row has `{TBD}`
      CIs.  See `docs/table3_executable.md` §0.1 for the item mapping.

## 5. What worker 050 did **not** do

- Did not push or upload anything; did not contact Zenodo/GitHub.
- Did not invent a repository URL or DOI.
- Did not change the manuscript, figures, tables, or the 00–03 ledgers.
- Did not change the algorithm or any saved result on `L-wdl`; the fixture was
  read from `dataset_v3` and written only to `/tmp` on the server, then copied
  here.
- Left the real-data fixture and all expected outputs in the package so the
  release can be checked without server access.
