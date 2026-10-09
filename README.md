# Stochastic Kinetic Screening (SKS)

A source-code implementation of **stochastic kinetic screening for catalytic reaction networks**.

The workflow samples unknown transition-state (TS) contributions over a prescribed energy range, evaluates an ensemble of CATKINAS microkinetic models, ranks kinetically influential elementary steps using degree-of-rate-control (DRC) occurrence statistics, and tests progressively larger fractions of explicitly refined reaction barriers.

> [!IMPORTANT]
> This repository contains the original research scripts rather than a fully portable software package. Several paths, scheduler identifiers, sampling parameters, and thresholds are hard-coded and must be checked before use.

---

## Method Overview

For elementary step $i$, the provisional forward activation free-energy barrier is assigned as

$$
\Delta G^\ddagger_{i,f}
=
\max(\Delta G_i,0)+E_{x,i},
$$

where $E_{x,i}$ is a non-negative sampled TS contribution.

The workflow contains two stages:

1. **Stochastic barrier sampling (`sks.py`)**
   - Generates Latin hypercube samples for $E_{x,i}$.
   - Creates multiple CATKINAS input files.
   - Copies each sampled input to `INCAR.m`.
   - Submits one MATLAB/CATKINAS calculation at a time through PBS/Torque.
   - Produces an ensemble of microkinetic output logs.

2. **Statistical screening and TS-refinement tests (`N/N_test.py`)**
   - Reads the stochastic CATKINAS logs.
   - Extracts $|X_{\mathrm{DRC},i}|$ for every elementary step.
   - Counts how frequently each step exceeds the selected DRC threshold.
   - Ranks elementary steps by occurrence count.
   - Generates input files retaining increasing fractions of the highest-ranked reference barriers.
   - Submits the corresponding CATKINAS calculations for pathway and rate-convergence tests.

---

## Repository Structure

The original scripts assume the following relative directory layout:

```text
.
├── sks.py                    # Stochastic sampling and ensemble job submission
├── CATKINAS.p                # CATKINAS MATLAB P-code*
├── INCAR1.m                  # Fully parameterized reference input
├── INCAR.m                   # Active CATKINAS input, overwritten during the workflow
├── kinetic.script            # PBS/Torque job script
├── INCAR_1.m                 # Generated stochastic inputs
├── INCAR_2.m
├── ...
├── INCAR_100.m
├── results_1/                # Ensemble CATKINAS output: log1, log2, ...
└── N/
    ├── N_test.py             # DRC statistics and progressive TS refinement
    ├── CATKINAS.p            # CATKINAS MATLAB P-code*
    ├── INCAR1.m              # Reference input used in the refinement tests
    ├── INCAR.m               # Active refinement input
    ├── kinetic.script        # PBS/Torque job script for refinement tests
    ├── INCAR_0.1.m           # Generated refinement inputs
    ├── INCAR_0.2.m
    ├── ...
    ├── INCAR_1.0.m
    └── results_1/            # CATKINAS outputs from the refinement tests
```

\* Include or redistribute `CATKINAS.p` only when you have permission to do so. Otherwise, instruct users to provide their own authorized CATKINAS installation.

### Why the directory layout matters

The original source uses relative paths. In particular:

- `sks.py` must be run from the repository root.
- `N_test.py` must be run from the `N/` directory.
- `N_test.py` reads the stochastic inputs through paths such as `../INCAR_1.m`.
- The stochastic CATKINAS logs are expected in the root-level `results_1/` directory.
- Each `kinetic.script`, `CATKINAS.p`, and active `INCAR.m` must be located in the same CATKINAS working directory.
- CATKINAS creates `results_1/` at the same level as `CATKINAS.p`.

---

## Requirements

### Python

Python 3.8 or later is recommended.

Install the required Python packages with:

```bash
pip install numpy scipy pyDOE2
```

The scripts also use the following Python standard-library modules:

- `os`
- `re`
- `sys`
- `time`
- `shutil`
- `threading`
- `subprocess`

### External software

The workflow additionally requires:

- MATLAB
- CATKINAS
- A PBS/Torque-like scheduler providing:
  - `qsub`
  - `qstat`
  - `qdel`
- A cluster-specific `kinetic.script` that starts MATLAB and runs CATKINAS in the current working directory

The exact MATLAB command, module-loading commands, CPU settings, wall time, queue, and account information must be adapted to the local cluster.

---

## CATKINAS Input Format

The scripts assume that the elementary-reaction block in `INCAR1.m` is located between:

```text
% examples :
```

and:

```text
Q0 = [1 1 1];
```

Each reaction should appear on one line, for example:

```text
CH2CO_in+#1<->CH2CO_H#1        [ 0.7500   -0.4000 ]
```

The two bracketed values are interpreted as:

```text
[ forward_activation_barrier   reaction_free_energy ]
```

For a stochastic realization, the forward barrier is replaced according to:

```python
Ea = max(deltaG, 0.0) + Ex
```

while the reaction free energy remains unchanged.

---

## 1. Stochastic Barrier Sampling

From the repository root, run:

```bash
python sks.py
```

The current source-code defaults are:

```python
N_REACTIONS = 81
n_samples = 100
X_max = 1.5
random_state = 43
```

`pyDOE2.lhs` is used with the `maximin` criterion to generate the sampling matrix.

For every realization, the script:

1. creates `INCAR_<sample_id>.m`;
2. copies it to the active input `INCAR.m`;
3. submits `kinetic.script`;
4. waits for the current CATKINAS calculation to finish;
5. proceeds to the next realization.

The calculations are intentionally submitted sequentially because all jobs use the same active file, `INCAR.m`.

---

## 2. DRC Statistics and Step Ranking

After all stochastic calculations finish, enter the refinement directory:

```bash
cd N
python N_test.py
```

`N_test.py` reads the stochastic CATKINAS output logs from the root-level ensemble-results directory and associates `log<n>` with `../INCAR_<n>.m`.

For every elementary step, the script records:

- DRC values across stochastic realizations
- sampled TS contributions
- mean DRC
- DRC standard deviation
- Pearson correlation between the sampled barrier and DRC
- linear-response slope
- the number of realizations in which the DRC magnitude exceeds the threshold

### Current DRC occurrence threshold

The provided source currently uses:

```python
abs(X_DRC) > 1.0e-4
```

This threshold is a study-specific parameter rather than a universal SKS constant. Modify it in `get_tof_infulence()` when another value is required.

### Ranking rule

Although several statistics are calculated, the final elementary-step ranking is based on the occurrence count:

```python
corr_list.sort(key=lambda x: x[4], reverse=True)
```

The analysis outputs include:

```text
corr_list.txt
drc_number.txt
```

---

## 3. Progressive TS Refinement

`N_test.py` tests the following retained-barrier fractions:

```python
threshold_list = [
    0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 1.0
]
```

For a fraction $N$:

- the top-ranked $N$ fraction of elementary steps retain their reference barriers from `N/INCAR1.m`;
- all remaining steps receive provisional barriers according to

$$
\Delta G^\ddagger_{i,f}
=
\max(\Delta G_i,0)+0.75\ \mathrm{eV}.
$$

The generated files are:

```text
INCAR_0.1.m
INCAR_0.2.m
...
INCAR_1.0.m
```

Each generated file is copied to `N/INCAR.m` and submitted through `N/kinetic.script`.

The corresponding CATKINAS results are written to:

```text
N/results_1/
```

These calculations can then be used to examine convergence of:

- the favored reaction pathway;
- the catalytic reaction rate.

---

## Required Manual Configuration

Before running the scripts on a new system, inspect and modify the following settings.

### 1. Stochastic-log directory in `N/N_test.py`

The source contains a project-specific absolute path:

```python
cwd = "/home/hxyang/.../result_1"
```

Replace it with the absolute path to the root-level stochastic output directory, for example:

```python
cwd = "/path/to/repository/results_1"
```

The directory must contain:

```text
log1
log2
...
```

### 2. `result_1` versus `results_1`

CATKINAS in the present setup creates:

```text
results_1/
```

Some original source-code lines may still use:

```text
result_1/
```

Check both scripts and replace the directory name consistently with the actual CATKINAS output directory before running.

### 3. Scheduler job identifier

The queue-monitoring commands currently search for:

```bash
grep curve_CO
```

Replace `curve_CO` with the job name defined in your own `kinetic.script`.

### 4. Number of reactions

In `sks.py`:

```python
N_REACTIONS = 81
```

This value must equal the number of elementary-reaction lines in the sampled reaction block.

### 5. Sampling parameters

Check the following source-code values:

```python
n_samples = 100
X_max = 1.5
random_state = 43
```

The loop:

```python
for sample_id in range(100):
```

must remain consistent with the selected number of samples.

### 6. DRC threshold

In `N/N_test.py`, check:

```python
if abs(value) > 0.0001:
```

Set this to the threshold used in the corresponding SKS analysis.

### 7. CATKINAS and MATLAB job script

Both root-level and `N/` copies of `kinetic.script` must be adapted to the local:

- MATLAB executable or module;
- queue name;
- project/account;
- requested CPU resources;
- wall time;
- CATKINAS invocation.

### 8. Queue monitoring and cleanup

The original scripts contain cluster-specific `qstat`, `qdel`, and failure-cleanup commands. Review these commands carefully before execution. In particular, confirm that they cannot cancel unrelated jobs or remove unrelated CATKINAS result files.

---

## Output Interpretation

A high occurrence count means that an elementary step repeatedly exhibits a DRC magnitude above the prescribed threshold across many stochastic barrier assignments.

Such a step is considered kinetically influential over a broader region of the sampled barrier space and is prioritized for explicit TS refinement.

This statistical interpretation differs from identifying a rate-controlling step from a single deterministic microkinetic model.

---

## Reproducibility Notes

The current Latin hypercube sampling uses a fixed random seed:

```python
random_state = 43
```

For exact reproducibility, retain:

- the original Python and dependency versions;
- all generated `INCAR_<n>.m` files;
- the corresponding `results_1/log<n>` files;
- the final ranking files;
- the refinement input and output files.

Saving only the random seed may not guarantee identical samples across different package versions.

---

## Methodological Summary

The SKS workflow combines:

1. thermodynamically consistent stochastic barrier sampling;
2. ensemble microkinetic simulations;
3. DRC occurrence statistics;
4. ranking of kinetically influential elementary steps; and
5. progressive DFT-level TS refinement.

The objective is not to select one arbitrary provisional kinetic model, but to identify the elementary steps that remain kinetically influential across an ensemble of possible unknown TS contributions.

---

## Citation

If this repository is used in academic work, please cite the associated SKS publication after its bibliographic information becomes available:

```text
Adaptive Stochastic Kinetic Screening for Identifying Catalytic Reaction Pathways
```

---

## License

No license is currently specified.

Before public release, confirm:

1. which license should apply to `sks.py` and `N_test.py`; and
2. whether `CATKINAS.p` may legally be redistributed in the repository.

A third-party or separately licensed CATKINAS file should not automatically be covered by the license applied to the Python scripts.

---

## Contact

For questions about the methodology or implementation, please open a GitHub issue or contact the repository authors.
