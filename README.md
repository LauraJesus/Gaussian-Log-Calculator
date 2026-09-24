# Gaussian Log Calculator

Extracts ESP and NBO atomic charges from many Gaussian output (`.log`) files, averages them across configurations, and closes the solute charge sum onto a target value by distributing the residue among chemically equivalent groups.

## Context
Built as part of my computational chemistry research. The two scripts cover the two halves of one step of an **iterative polarization cycle**: a Monte Carlo simulation produces dozens or hundreds of solute configurations, each one submitted to Gaussian; the charges of every configuration have to be averaged, and the resulting average has to be turned back into a charge set whose sum is physically correct before it can feed the next cycle. Doing that by hand across hundreds of logs is unfeasible.

```
Gaussian .log files ──► scriptLog.py ──► *_cargas_medias.txt ──┐
                                                               ├──► ajuste/ajustev2.py ──► DICE file for the next cycle
initial DICE file (ljname) ────────────────────────────────────┘
```

## Requirements
- Python 3.9+
- No external dependencies — standard library only (`os`, `re`, `sys`, `math`, `argparse`, `pathlib`, `decimal`, `dataclasses`)

---

# 1. `scriptLog.py` — charge extraction and averaging

## What it does
- Reads every `.log`/`.txt` file in a folder, sorted **numerically** by filename (`log2`, `log10`, `log20` in the right order, not alphabetically)
- Extracts ESP charges from the `ESP charges:` block and NBO charges from the `Summary of Natural Population Analysis:` block
- Flags files that are incomplete or failed: crash/timeout, Gaussian error termination, out of memory, SCF convergence failure
- Optionally discards individual charges that violate a **sign restriction** declared by the user — a configuration in which, say, a chlorine comes out positive is dropped for that atom only, not for the whole file
- Computes N (how many configurations actually entered the average), mean, and standard deviation for each solute atom
- Reports the sum of the solute's mean charges at the end

## Usage
```bash
python scriptLog.py <folder> --solute <number_of_solute_atoms> --output <output_folder>
```

With no arguments the script asks for the folder and the solute size interactively:
```bash
python scriptLog.py
```

### Options
- `--solute N` — number of atoms in the solute (the first N atoms of each log)
- `--output DIR` — output folder (default: `results`)
- `--sem-outliers` — skip the outlier questions entirely (non-interactive runs)
- `--decimais N` — decimal places in the means TXT (default: 6)

### Outlier removal
If not disabled, the script asks — separately for ESP and for NBO — whether to apply sign restrictions. Each rule targets either an element symbol (`Pt`) or an atom selection (`2,3`, `4-7`), and sets a bound (`> value` or `< value`). Charges violating the rule are discarded per configuration, and the script reports how many were dropped for each atom. An atom that loses every configuration ends up with N = 0 and a mean of `nan`, flagged in the output.

## Output
In the output folder:

- **`esp_cargas_table.html`** / **`nbo_cargas_table.html`** — full table for visual inspection: atom number, symbol, N, mean, standard deviation, and one column per processed file, with failed files and discarded charges highlighted
- **`esp_cargas_medias.txt`** / **`nbo_cargas_medias.txt`** — one line per solute atom, `symbol mean`, in atom order. This is the file consumed by the adjustment step

---

# 2. `ajuste/ajustev2.py` — charge residue adjustment

The averages coming out of step 1 do not sum to the solute's formal charge: averaging over configurations leaves a residue. This script closes that sum while keeping chemically equivalent atoms identical.

## Method
Equivalent atoms receive the same adjusted charge. For each group the script computes the width of the interval between its initial charge and its mean charge, and distributes the residue **proportionally to those widths** — a group whose charge barely moved between cycles absorbs little of the correction, while a group that moved a lot absorbs more.

```
residue   = target − Σ (n_g · q_mean_g)
f         = residue / Σ (n_g · width_g)
q_final_g = q_mean_g + f · width_g
```

Equivalence groups are **declared by the user**, not derived from automatic bond perception. This follows the RESP practice (Bayly, Cieplak, Cornell & Kollman, *J. Phys. Chem.* **1993**, 97, 10269), where equivalence restraints are imposed by the chemist. An earlier version inferred connectivity geometrically, from covalent radii times a 1.3 tolerance factor; that factor is an implementation choice with no published basis, so it was dropped in favor of an explicit declaration.

The same groups file serves both the ESP and the NBO tracks — grouping depends only on the solute topology, never on the charges. Keeping the partition identical across both models is what makes the difference between them reflect the charge model alone.

## Inputs
1. **Initial DICE file (`ljname`)** — the one used in the cycle's first simulation (e.g. `s0Pt001DFT.txt`). Initial charges are read from the charge column of the first molecule block (the solute, by DICE convention). Everything else in the file is kept verbatim for the output.
2. **Means TXT** produced by `scriptLog.py` (`esp_cargas_medias.txt` / `nbo_cargas_medias.txt`). The charge model (ESP or NBO) is detected from the filename.
3. **Groups file** declaring chemical equivalence, one `label : selection` per line. `#` starts a comment. Numbering follows the solute block of the DICE file:

```
# c,t-[Pt(en)Cl2(OH)2] (Pt001)
Pt          : 1
Cl          : 2,3
O_hidroxo   : 4,6
H_hidroxo   : 5,7
N_amino     : 8,10
H_amino     : 9,11-13
C_metileno  : 14,15
H_metileno  : 16-19
```

The partition must be total and exclusive — every solute atom in exactly one group. A missing or repeated atom raises an error, so a typo stops the run instead of silently producing a wrong grouping.

## Usage
```bash
python ajuste/ajustev2.py --iniciais s0Pt001DFT.txt \
                          --medias results/esp_cargas_medias.txt \
                          --grupos ajuste/grupos_Pt001.txt \
                          --atomos todos
```

Any argument left out is asked for interactively.

### Options
- `--iniciais FILE` — initial DICE (`ljname`) file
- `--medias FILE` — means TXT from `scriptLog.py`
- `--grupos FILE` — equivalence groups file
- `--sem-grupos` — impose no equivalence, each atom its own group; diagnostic use only
- `--atomos SEL` — which atoms take part in the adjustment: `1-19`, `1,4,9`, `1-13,15,17-19`, or `todos`. A selection that splits an equivalence group is rejected
- `--alvo Q` — desired total solute charge (default: 0)
- `--criterio-grupo {media,menor_largura}` — how each group's `(q_initial, q_mean)` pair is defined: `media` takes the arithmetic mean over the group's atoms; `menor_largura` takes the single atom with the smallest individual `|q_mean − q_initial|`, which avoids sign cancellation inside the group. Asked interactively if omitted. Neither criterion has a published precedent
- `--decimais N` — decimal places in the table and in the written files (default: 6)

## Output
The adjustment table is always printed to the terminal: a per-group block (n, initial charge, mean charge, width, displacement, final charge) followed by a per-atom block, and the sum over all atoms as a check. Two files can be written, each asked for after the table is shown, or passed directly for automation:

- **`--saida FILE`** — a copy of the input DICE file with the solute charges replaced by the adjusted ones. Only the charge field of each solute line changes; labels, atomic numbers, coordinates, Lennard-Jones parameters, original spacing, title, solvent block and `$end` stay byte-for-byte identical, so the file goes straight into the next cycle
- **`--saida-txt FILE`** — two columns, `q_mean q_final`, one line per solute atom, for plotting or record keeping

After writing the DICE file the script sums the charges **as rounded on disk** and warns if rounding left a residue above 1e-5 — raising `--decimais` fixes it.

---

## Repository layout
```
scriptLog.py              charge extraction and averaging
ajuste/
  ajustev2.py             residue adjustment
  grupos_Pt001.txt        example groups file
```

### Folders you need to provide
Only the code is versioned. The data a real run consumes and produces is not — Gaussian logs alone run into hundreds of files per cycle — so these folders have to exist locally:

| Folder | Contents | Used by |
|---|---|---|
| `cargaschelpg/` | Gaussian `.log` (and matching `.gjf`) files of the ESP/CHELPG calculations, one pair per Monte Carlo configuration | input to `scriptLog.py` |
| `cargasnbo/` | the same, for the NBO calculations | input to `scriptLog.py` |
| `valores/` | the initial DICE (`ljname`) files of each cycle, e.g. `s0Pt001DFT.txt`, `s0Pt001D_NBO.txt` | input to `ajustev2.py` (`--iniciais`) |
| `results/` | tables and means written by `scriptLog.py` | output, then input to `ajustev2.py` (`--medias`) |

A DICE file in `valores/` looks like this — the combination rule, the number of molecule types, the title line, and then one line per atom (number, atomic number, x, y, z, charge, epsilon, sigma), the solute being the first block:

```
*
2
19 JFL (2009) Pt001 M062X CHELPG
1         78      -0.189864    0.000289   -0.034045      0.412134   7.0100   2.559
2         17      -0.239781    0.594002    2.230323     -0.325673   1.7913   3.460
...
```

`ajustev2.py` reads the charge column of that first block and rewrites only that column in its output.

## Limitations / next steps
- Assumes a specific Gaussian output format (the `ESP charges:` and `Summary of Natural Population Analysis:` blocks) — it does not generalize to every calculation type
- No automated tests
- `scriptLog.py` does not assemble a complete `ljname`: coordinates and Lennard-Jones parameters are not present in Gaussian logs, which is why the adjustment step rewrites an existing DICE file instead of building one
- The `menor_largura` criterion and the width-proportional distribution itself are choices made for this project, without published precedent — they are documented here so the results can be read with that in mind
