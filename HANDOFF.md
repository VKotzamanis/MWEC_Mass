# HANDOFF — WaveConditions wave-climate plots

**Branch:** `claude/wave-station-plots-dozc0n` (3 commits ahead of the uploaded baseline)
**Scope:** `WaveConditions/` only. Nothing outside that folder was touched.
**Status:** Complete and pushed, but **never executed** — see [Critical caveat](#critical-caveat).

---

## Goal

The user (Vasileios) needed publication-quality plots of the **wave station data by
itself** — the WIS wave climate, with no wave-energy-converter in the picture. Their
words: *"The plots I need us to improve/create are simple and they have to do ONLY
with the wave station data themselves."*

Two concrete asks arrived over the session:

1. **W1–W4** — general station-only figures (occurrence, spectrum, energy, cross-station).
2. **W5** — reproduce a slide of theirs (three stations side by side, occurrence
   histograms over period on top and over `Hs` below, Times New Roman, bold location
   header with water depth) but on **`Te` instead of `Tp`**.
3. **W6** — a 3-D `surf` per station, **X = Te, Y = Hs, Z = probability-weighted
   available wave energy**, "as smooth as possible", using the spectral formula.

Quality bar was stated explicitly: *"The plots would need to be a masterpiece quality."*
Target is **MATLAB R2021b**.

---

## Where the WIS data actually comes from (read this first)

`WaveConditions/run_buildClimateGrid.m` is the driver that parses WIS files and writes
`<station>_climate_grid.mat`.

**Its worker functions are NOT in the repo.** Line 42 adds a folder
`WaveConditions/MATLAB Script for Plotting/` that was never committed. That folder holds
`parseWisOneline`, `parseWamSpectra`, `buildClimateGrid`, `writeClimateGridMat` and
`plotClimateGrid`. So:

- `run_buildClimateGrid.m` **cannot run** as checked out.
- `plotClimateGrid` — the pre-existing station-only plotting — is among the missing files.
  That is almost certainly why this work was requested.
- Only the four **outputs** are committed.

Do not waste time looking for those functions; they are not in git history either.

### The four climate grids

`WaveConditions/WIS_Output_WAM/<station>_climate_grid.mat`, top-level var `climateGrid`,
schema `STABILITY_HANDOFF_PLAN v2.1.0` (empirical-only: every occupied cell stores a
measured mean spectrum; no JONSWAP/Bretschneider closure anywhere).

| station | region (stored) | depth | Hs bins | Te bins | N records |
|---|---|---|---|---|---|
| ST63044 | `NorthAtlantic` | 74 m | 6 (0–3.0) | 13 (0–13) | 8685/8772 |
| ST73135 | `SouthPass` | 28 m | 6 (0–3.0) | 13 (0–13) | 8696/8784 |
| ST73329 | `Mexico` | 427 m | 8 (0–4.0) | 12 (0–12) | 8751/8784 |
| ST84040 | `Pacific` | 50 m | 4 (0–2.0) | 16 (0–16) | 8750/8772 |

**Native intervals: dHs = 0.50 m, dTe = 1.00 s.** Each cell carries a 115-point spectrum
over ω = 0.30…6.00 rad/s at dω = 0.05.

Array orientation (verified against the files, easy to get backwards):
- `probability_grid` is **nHs × nTe**; `S_omega_grid` is **nHs × nTe × nOmega**.
- `Hs_centers` is nHs×1, `Te_centers` is nTe×1, `omega` is nOmega×1.
- `source_grid` and `fit_reason_grid` are MATLAB **string arrays** stored as opaque MCOS
  objects in the v7.3 file — unreadable from h5py, and not needed. Use
  `probability_grid > 0` plus a finiteness check on the spectra instead.
- `iec62600_101_summary.P_wave_kWm` is **NaN** and `.epsilon` is **0** in all four files.
  Both are recomputed in code; do not trust the stored values.

---

## What was built

Two new files, both in `WaveConditions/`:

### `MWEC_WaveClimate_Plots.m` (915 lines, classdef, all static)

| line | member | what |
|---|---|---|
| 33 | `style` | Wong palette from `MWEC_Tuning_Plots.style`, Times New Roman swapped in |
| 47 | `serif_font` | Times New Roman → metric-compatible fallback, cached in a `persistent` |
| 72 | `derive(cg, opts)` | **the core** — Sew, J, per-cell E, P_wave, IEC moments, 90% band |
| 161 | `group_velocity` | deep vs finite-depth `c_g` |
| 186 | `wave_number` | solves ω² = gk·tanh(kh), Guo (2002) start + Newton |
| 245 | `fig_w1_scatter` | (Hs,Te) occurrence table + marginals, log colour |
| 331 | `fig_w2_spectrum` | Sew stacked by contributing sea state, + J(ω) + cumulative |
| 416 | `fig_w3_energy` | per-cell share of resource + energy-concentration curve |
| 483 | `fig_w4_stations` | cross-station spectra + one bar panel per metric |
| 556 | `fig_w5_histograms` | **the slide figure**, station per column, hand-positioned axes |
| 655 | `fig_w6_energy_surface` | **the 3-D surface** |
| 736 | `energy_density` | smoothing kernel that builds W6's Z |
| 768 | `reflected_gaussian` | Gaussian folded about the origin |
| 795 | `label_cells` | in-cell numerals, text colour flipped by colormap position |
| 820 | `cividis` | embedded 256-entry table |

### `run_wave_climate_plots.m` (153 lines)

Driver. Globs `WIS_Output_WAM/*_climate_grid.mat`, renders W1/W2/W3/W6 per station then
W5/W4 across stations, prints a resource table, saves `WaveClimate_summary.mat`.

```matlab
station_labels = struct('ST63044','North Atlantic, NH', 'ST84040','Santa Barbara, CA');
run_wave_climate_plots
```

CONFIG block at the top (all optional): `climate_dir`, `fig_dir`, `stations`, `formats`,
`sigma_factor`, `depth_model`, `station_labels`.

---

## The two pieces of physics that matter

### 1. Energy flux — now finite-depth

Per-cell omnidirectional flux, from the cell's **stored empirical spectrum**, never the
parametric `ρg²Hs²Te/(64π)`:

```
J_ij = ρ g ∫ c_g(ω) S_ij(ω) dω          [W/m]
c_g  = ½(1 + 2kh/sinh 2kh)·ω/k,   ω² = g k tanh(kh)
```

The code originally used deep water `c_g = g/(2ω)`, matching
`MWEC_Tuning_Kernels.build_S_ew:97`. **That is wrong for two of these stations.** Deep
water needs kh > π; kh at the energy period is **2.93 at South Pass (28 m)** and **3.87 at
Pacific (50 m)**. Switching to the full dispersion relation raises the resource:

| station | deep | finite | Δ |
|---|---|---|---|
| ST73135 South Pass (28 m) | 2.786 | **2.996** | +7.0% |
| ST84040 Pacific (50 m) | 1.599 | **1.706** | +6.3% |
| ST63044 N. Atlantic (74 m) | 3.156 | **3.210** | +1.7% |
| ST73329 Mexico (427 m) | 5.917 | 5.917 | +0.0% |

`opts.depth_model = 'deep'` restores the old form. The driver prints **both** numbers
every run so the assumption is never invisible. Spectra and IEC moments are unaffected —
only J, E and P_wave move.

> **Open item, deliberately not acted on:** `MWEC_Tuning_Kernels.build_S_ew` still uses
> the deep-water form and therefore understates `F_ew` by ~7% at the two shallow
> stations. The W figures now disagree with the T figures by those percentages. Changing
> it would shift tuning results, so it was left to the user. Raise it again if tuning
> numbers are ever compared against W-figure numbers.

### 2. W6's Z axis — a density, not a per-cell number

Native grid is only 4–8 × 12–16 cells, so a raw `surf` is a staircase. `energy_density`
convolves the cell weights `p_ij·J_ij` with a Gaussian of one bin width, **reflected about
Hs = 0 and Te = 0** so no energy leaks to negative values, evaluated on 181×241.

```
Z(Hs,Te) = Σ_ij p_ij·J_ij · K_H(Hs; Hs_i)·K_T(Te; Te_j)     [kW/m per (m·s)]
∬ Z dHs dTe = Σ p_ij J_ij = P_wave                          [kW/m]
```

**Why a density:** a per-cell value is mesh-dependent — halve the bin width and it halves.
`Z = ∂²P_wave/∂Hs∂Te` is mesh-independent, so refining converges instead of shrinking. The
invariant is the **volume**, not the peak height. Convolution with a normalised kernel
preserves the integral exactly, which is why the W6 subtitle prints the volume next to
`Σp_ij·J_ij` as a live self-check — they agree to ~0.03%.

Z is *available* (incident) flux — no device, no capture width, no PTO — and a time-mean,
since `p_ij` is a record fraction. × 8766 h gives kWh/m/yr.

`opts.smooth = false` surfs the native cells instead (units then kW/m per cell).

---

## What worked

- **Verifying MATLAB physics in Python.** No MATLAB in the container, so `derive` was
  ported to numpy and run against all four real grids. Caught the orientation question
  immediately and confirmed: m0 budget matches the Hs² target to **0.00%** at every
  station, `Σ(p·J_cell)` reconciles with `trapz(J)`, surface volume matches to 0.03%, and
  Guo+3-Newton hits the exact dispersion root to **2.2e-16**.
- **Rendering matplotlib mock-ups before writing MATLAB.** Every figure was prototyped and
  looked at first. This caught four design faults that would otherwise have shipped —
  see below.
- **Liberation Serif for previews** — metric-compatible with Times New Roman, so preview
  layout is faithful.
- **Single-sourcing style.** `MWEC_WaveClimate_Plots.style()` calls
  `MWEC_Tuning_Plots.style()` and overrides only the font, so geometry stays shared and
  the T figures were never touched.
- **Static balance checking** with a small Python paren/block counter, since there is no
  `mlint`. Note it mis-parses MATLAB's `''` escape — one flagged "imbalance" was a false
  positive.

## What didn't work

- **`apt-get install octave`** — apt is offline in this container. No way to syntax-check
  MATLAB here. Don't retry.
- **Overlaying raw per-cell spectra on W2.** `S_ew` is a probability-weighted mean, so it
  is smaller than *any* individual cell by construction and vanished under the family of
  curves. Replaced with a **stacked decomposition** (top 5 energy cells + remainder) whose
  areas sum exactly to `S_ew`.
- **Dark-at-top colormap for W3.** Made near-zero cells indistinguishable from unoccupied
  ones. Fixed with light-at-zero, then superseded by cividis.
- **Floor contour projection under the W6 surface.** Tried, then removed: the surface's own
  near-zero skirt spans the whole (Te,Hs) plane and occludes anything below it from every
  useful viewing angle. Don't re-add it.
- **Mixing units on one bar axis in W4** (m, s, kW/m) — only the largest unit stayed
  legible. Split into one panel per metric.
- **33-anchor cividis + interpolation** — ~2% max error concentrated in the blue ramp,
  visible, and breaks colour agreement with matplotlib. Embedded the full 256-entry table.
- **KDE overlay at σ = 0.75·dT on W5's Te bars** — visibly undershot the mode and read like
  a bad fit. Dropped to σ = 0.50·dT.

## Known-fragile / traps

- **`label_cells` colour polarity.** cividis is dark-low/light-high. If anyone swaps in a
  dark-at-top colormap, the white-numeral rule at line ~795 must flip or numerals go white
  on yellow.
- **W5 axes are hand-positioned** (`axes(fig,'Position',...)`, not `tiledlayout`) so the
  columns align and the rotated `Hs` interval labels get room. **This is the most likely
  thing to need nudging on first real run** — constants live at
  `fig_w5_histograms`, `L/R/gap/y_te/y_hs/rowh`.
- **`style()` self-recursion.** It must call `MWEC_Tuning_Plots.style()`. A blanket
  find-and-replace of `MWEC_Tuning_Plots.style()` → `MWEC_WaveClimate_Plots.style()` will
  make it call itself. This happened once and was caught.
- **R2021b compatibility.** Use `caxis` not `clim`; no `fontname()`; `subtitle` and
  `tiledlayout` are fine. Keep the `'tex'` interpreter — `'latex'` silently swaps in
  Computer Modern and breaks the Times match.
- **Station labels.** W5 headers default to the stored region. The user's slide reads
  Gulf Coast LA / North Atlantic NH / Santa Barbara CA at 103/74/50 m. Only two map here:
  74 m = ST63044, 50 m = ST84040. **Nothing in the repo is 103 m** (others are 28 m and
  427 m), so their Gulf Coast column came from a station not in `WIS_Output_WAM`. Set
  `station_labels` explicitly rather than guessing.

---

## Next steps

1. **Run it.** `run_wave_climate_plots` in `WaveConditions/`. This is the top priority —
   nothing has ever been executed. Expect to nudge W5's hand-positioned axes.
2. **Check the console self-checks** the driver prints per station: surface volume vs
   `Σp·J` (should be <0.1% apart) and finite vs deep-water P_wave. If either is off,
   something is wrong in `derive` or `energy_density`, not in the figures.
3. **Decide the `build_S_ew` question** — whether `MWEC_Tuning_Kernels` should also move to
   finite-depth `c_g`. Affects tuning results; user's call.
4. **Finer Te histogram** (if wanted): W5's bars are the grid's own 1.0 s bins. Per-record
   `Te` is discarded at grid-build time, so this needs a smaller `dbinTe` at the top of
   `run_buildClimateGrid.m` **and** the raw WIS files **and** the missing
   `MATLAB Script for Plotting/` parser folder. Ask the user for that folder first.
5. **Not requested, do not do unasked:** touching the T figures, opening a PR, or changing
   the tuning pipeline.

## Critical caveat

**None of this MATLAB has been executed.** There is no MATLAB and no Octave in the
container and apt is offline. What *is* verified:

- all physics, ported to numpy and checked against the four real grids (numbers above)
- paren/block balance, statically
- figure layouts, via matplotlib renderings of the same data

What is **not** verified: that the MATLAB parses, that handle-graphics calls behave as
intended, and that the layouts survive real font metrics. Treat the first run as a
debugging session, not a smoke test.
