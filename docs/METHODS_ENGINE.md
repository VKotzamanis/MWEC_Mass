# Methods and computational engine

This document describes the algorithms that the shipped MATLAB source actually executes. It is an
implementation reference for reviewers, analysts, and maintainers—not a proposal for a future model.
Where the code uses an approximation, fallback, or diagnostic calculation, that status is stated
explicitly.

The source remains authoritative for numerical tolerances and release-specific defaults. Runtime setup,
input fields, and output selection belong in `RUNTIME_GUIDE.md`; this document focuses on geometry,
physics, optimisation, and realisation.

## Contents

1. [Scope and computational sequence](#1-scope-and-computational-sequence)
2. [Coordinates, notation, and units](#2-coordinates-notation-and-units)
3. [MultiSurf ingestion and `MS2Parser`](#3-multisurf-ingestion-and-ms2parser)
4. [Production hydrostatics: contour tables and strip integration](#4-production-hydrostatics-contour-tables-and-strip-integration)
5. [Divergence theorem, Gauss–Legendre quadrature, and clipping](#5-divergence-theorem-gausslegendre-quadrature-and-clipping)
6. [Mass properties and flotation](#6-mass-properties-and-flotation)
7. [Stability, stiffness, hydrodynamics, and periods](#7-stability-stiffness-hydrodynamics-and-periods)
8. [Optimization formulation](#8-optimization-formulation)
9. [Realization models](#9-realization-models)
10. [Numerical safeguards and failure behavior](#10-numerical-safeguards-and-failure-behavior)
11. [Assumptions and limitations](#11-assumptions-and-limitations)
12. [Source map for the implemented method](#12-source-map-for-the-implemented-method)

## 1. Scope and computational sequence

The engine converts a MultiSurf hull definition and a cached hydrodynamic data set into an optimised
vertical mass distribution and, optionally, an as-built material realisation.

```text
MultiSurf .ms2 geometry
        │
        ├─> evaluatable parametric surfaces and cached iso-z boundaries
        │
        ├─> adaptive waterplane tables A_w(z), I_wp(z), P(z)
        │         └─> cumulative V_sub(z), CB_z(z), strip geometry
        │
hydrodynamic cache ─> interpolate A∞ and B̄ at vertical shift; transfer to candidate CG
        │
        └─> Stage 1 warm start ─> Stage 2 density optimisation
                                      │
                                      └─> preliminary | thin_shell | modular_precast
```

The production optimisation path is table based. It obtains displaced volume and buoyancy centre by
trapezoidal integration of cross-sectional area, not by repeatedly integrating clipped parametric
surfaces. Divergence-theorem/Gauss–Legendre routines remain important for initial global geometry,
validation, wetted area, and table-absent fallbacks; Section 5 distinguishes those roles precisely.

Principal entry points:

- [driver/run.m](../src/+mwecmass/+driver/run.m) — top-level pipeline.
- [driver/build_config.m](../src/+mwecmass/+driver/build_config.m) — normalised configuration.
- [optim/run.m](../src/+mwecmass/+optim/run.m) — Stage 1 and Stage 2.
- [realise/preliminary/run.m](../src/+mwecmass/+realise/+preliminary/run.m),
  [realise/thin_shell/run.m](../src/+mwecmass/+realise/+thin_shell/run.m), and
  [realise/modular_precast/run.m](../src/+mwecmass/+realise/+modular_precast/run.m) — realisation dispatch.

## 2. Coordinates, notation, and units

### 2.1 Frames and sign conventions

The geometry is defined in a body-fixed Cartesian frame `(x_b,y_b,z_b)`. The evaluation frame has its
undisturbed free surface at `z=0`; `z` is positive upward in both frames. The design variable
`vertical_shift`, denoted `v_s`, translates the body vertically:

\[
z = z_b + v_s, \qquad z_{wl,b}=-v_s.
\]

A positive `v_s` moves the hull upward. Consequently, the waterline moves downward in body coordinates
and physical draft decreases. With body-frame keel elevation `z_min`, the implementation reports

\[
d=\left|z_{min}+v_s\right|.
\]

The absolute value is safe inside the configured shift bounds, which keep the keel below the free
surface. It should not be interpreted as a general signed immersion formula outside those bounds.

The optimisation design vector is

\[
\mathbf{x}=[v_s,\rho_1,\ldots,\rho_N]^T,
\]

where strip 1 is the lowest horizontal strip and strip `N` is the highest. The dynamic three-degree-of-
freedom order is `[surge, heave, pitch]`; six-degree-of-freedom rigid-body matrices use
`[surge, sway, heave, roll, pitch, yaw]`.

### 2.2 Symbol table

| symbol | code name | meaning | SI unit |
|---|---|---|---|
| `v_s` | `vertical_shift` | body-to-world vertical translation | m |
| `z_wl,b` | `z_wl` | waterline elevation in body frame | m |
| `A_w` | `Aw` | waterplane area | m² |
| `I_wp,xx`, `I_wp,yy` | `I_wp_xx`, `I_wp_yy` | waterplane second moments about the body axes | m⁴ |
| `V_sub` | `V_sub` | submerged/displaced volume | m³ |
| `CB_z`, `CG_z` | `CB(3)`, `CG_total(3)` | buoyancy and mass-centre elevations | m |
| `BM`, `KM`, `GM` | same | longitudinal metacentric quantities for pitch | m |
| `M` | `mass_total` | total physical mass | kg |
| `M_b` | `mass_buoyant_force` | displaced water mass, despite the historical field name | kg |
| `A∞` | `A_full` | infinite-frequency added-mass matrix | mixed generalized units |
| `B̄` | `B_full` | band-averaged radiation-damping matrix | mixed generalized units |
| `K` | `K_total` | hydrostatic stiffness; PTO contribution is zero | mixed generalized units |
| `T_h`, `T_p` | `periods.heave`, `.pitch` | coupled-mode natural periods | s |
| `rho_w` | `RHO_WATER` | water density | kg/m³ |
| `g` | `G` | gravitational acceleration | m/s² |

`mass_buoyant_force` is a legacy name: the stored quantity is `rho_w V_sub` in kilograms, not force in
newtons. Force or stiffness formulas multiply it—or the equilibrated physical mass—by `g` explicitly.

Source: [hydrostatics/properties_3d.m](../src/+mwecmass/+hydrostatics/properties_3d.m).

## 3. MultiSurf ingestion and `MS2Parser`

### 3.1 File and record syntax consumed by the parser

The parser reads header metadata before `BeginModel;`, then entity records through `EndModel;`.
Recognized header fields are:

- `Units:` — stored as metadata;
- `Extents:` — six header bounds `[xmin ymin zmin xmax ymax zmax]`;
- `Symmetry:` — any listed `x` and `y` reflection planes.

An entity record has the general shape

```text
EntityType EntityName <entity metadata ...> / <type-specific references and values ...> ;
```

Token 1 selects the entity parser, token 2 is the entity name, and token 4 supplies the numeric
visibility value used by the implemented entity readers. `/` separates MultiSurf metadata from the
operands used by this implementation. Braced lists name control points, beads, component curves, or
loft sections. Indented lines and lines beginning with `{` or `A:` are appended to the current record
until its terminating semicolon. `Attribute:` lines within the model are ignored.

This is a supported subset, not a general MultiSurf grammar. `RealList`, `Variable`, and `Pathname`
records are deliberately skipped. An unknown entity type has no switch branch and is silently ignored;
a recognized record that throws while parsing emits a `ParseError` warning. A required entity that was
ignored will normally fail later during dependency resolution or evaluation.

> **Units warning:** `Units:` is recorded as parser metadata but coordinates are not converted. The field
> can be displayed by `geometry_summary`; it is not printed as part of the normal parser banner. The
> numerical engine assumes metre-based geometry. A deck declared in feet must be converted before use.

Source: [geometry/MS2Parser.m](../src/+mwecmass/+geometry/MS2Parser.m).

### 3.2 Render-supported physical surfaces

Only the following five types become evaluatable physical surfaces (`is_surface=true`). All use
normalized parameters `u,v in [0,1]`.

| surface | implemented definition | important behavior |
|---|---|---|
| `RuledSurf` | `S(u,v)=(1-v)C1(u)+vC2(u)` | analytical `S_u,S_v`; both boundaries share the same `u` parameterization |
| `RevSurf` | profile rotated around a line over a stored angle interval | Rodrigues-style axial/radial construction; axis cached; degrees converted to radians |
| `BLoftSurf` | B-spline in `v` through section curves evaluated at common `u` | open, clamped, uniformly generated knot vector; section parameterizations must correspond |
| `DevSurf` | `S(u,v)=(1-v)Snake(u)+vCurve(u)` | implemented as a direct ruled blend, not an independent developability solver |
| `MirrSurf` | reflected source surface | only `X=0` and `Y=0` surface mirrors are implemented |

For a revolution surface, let `O` be the axis origin, `a` its unit direction, and `P(u)` the profile.
With `P-O=z_a a+r`, and `e_t=a×e_r`, the evaluator uses

\[
S(u,v)=O+z_a a+\lVert r\rVert[\cos\phi(v)e_r+\sin\phi(v)e_t].
\]

At the axis (`||r||<10^-12`) it returns the profile point and a zero `v` derivative. `BLoftSurf`
evaluates the section curves first, then uses those points as control points of the loft-direction
B-spline. Because the knots are synthesized by this parser rather than imported from a complete native
MultiSurf representation, loft geometry should be independently checked when exact CAD reproduction is
required.

### 3.3 Supporting point, curve, and snake entities

The evaluator supports these non-surface dependencies:

| category | entities |
|---|---|
| points / point-on-parent | `FramePoint`, `MirrPoint`, `AbsBead`, `AbsRing` |
| curves | `BCurve`, `Conic`, `CopyCurve`, `Line`, `BSubCurve`, `Arc`, `PolyCurve2`, `ProjCurve` |
| snakes | `EdgeSnake`, `BSubSnake` |

`BCurve` uses a self-contained Cox–de Boor evaluator with an open clamped knot vector. `PolyCurve2`
allocates an equal normalized parameter interval to each component curve; it is not arc-length
parameterized. `BSubCurve` and `BSubSnake` map `[0,1]` between the first and last bead parameters and do
not construct a new spline from intermediate beads. `ProjCurve` zeroes the coordinate normal to the
named coordinate plane. `EdgeSnake` uses this surface-edge convention:

| edge | surface parameterization |
|---:|---|
| 1 | `v=0`, increasing `u` |
| 2 | `u=1`, increasing `v` |
| 3 | `v=1`, increasing `u` |
| 4 | `u=0`, increasing `v` |

Some native metadata is parsed but does not alter the implemented evaluator. For example, `Arc` stores
`arc_type`, `Conic` stores `conic_type`, and subcurve entities store a degree, while the live formulas
are selected by the entity class itself. Such records are supported only to the extent described above.

Dependencies are resolved recursively by name. `get_required_entities` performs a depth-first traversal
and returns parents before children. The point cache memoizes `FramePoint` and `MirrPoint` evaluations
only; a separate cache stores resolved revolution axes. Other point-on-parent, curve, snake, and surface
evaluations are not generally memoized. Mirror-chain resolution has a 20-level cycle guard.

### 3.4 File symmetry and mirror synthesis

If the model contains no explicit visible `MirrSurf` and its header declares symmetry, the parser creates
synthetic mirror surfaces. It first mirrors every original visible surface in `X` when requested. It then
mirrors the current visible set in `Y`, which also produces the combined `X+Y` reflection. If any explicit
visible `MirrSurf` exists, automatic header-symmetry synthesis is disabled for the entire model.

Reflections are not re-integrated unnecessarily. Surface classification resolves each mirror to its
ultimate source; volume and second moments are copied because they are reflection invariant, while the
appropriate centroid component changes sign.

### 3.5 Sampling, profiles, and horizontal contours

The parser remains parametric, but several production tables use sampled helper representations:

- actual hull `z` limits are estimated on a `50×50` parameter grid per visible surface rather than from
  header/control-point extents;
- the 2-D midplane profile samples surface boundaries at 200 parameter values and keeps points within
  `|y|<0.05 m`; it then mirrors `x` to form a bilateral profile;
- the production boundary cache stores 100 `u` samples for each source surface;
- contour points are angularly ordered about their arithmetic mean before polygon formulas are applied.

Horizontal iso-curves are dispatched by physical surface type:

- `RuledSurf` and `DevSurf`: solve the linear blend directly at each cached `u`,
  `v=(z_wl-z1)/(z2-z1)`. A coincident horizontal pair uses `v=0.5` only when it lies on the plane.
- `RevSurf`: find the first bracketed profile crossing in `u` by bisection, then analytically sample the
  stored revolution-angle arc.
- `BLoftSurf`: detect an azimuthal loft and sweep its loft parameter at each profile crossing; otherwise,
  at each `u`, solve the first bracketed loft-direction B-spline crossing by bisection.
- `MirrSurf`: reflect the already computed ultimate-source contour.

Non-monotone revolution profiles and axial B-spline lofts may have multiple crossings; the helper returns
the first. The B-spline helper warns when its span endpoints reveal multiple sign changes. Multiple
disconnected waterplane loops are not represented as separate polygons: all returned contour points are
combined and angularly ordered. These are material limitations for non-convex or multiply connected hulls.

C1 coverage is intentionally narrower than the parser's full implemented subset. The shipped C1 deck
contains literal `RevSurf` and `RuledSurf` surfaces plus supporting `FramePoint`, `MirrPoint`, `BCurve`,
`AbsBead`, `BSubCurve`, `Line`, `Arc`, `PolyCurve2`, `EdgeSnake`, and `ProjCurve` entities. Its header
causes synthetic `X` and `Y` mirrors. C1 does **not** exercise literal `BLoftSurf`, `DevSurf`, or
`MirrSurf` records, nor `Conic`, `CopyCurve`, `AbsRing`, or `BSubSnake` dependencies. Those types are
implemented, but C1 results are not evidence of their geometric fidelity.

Sources:
[geometry/precompute_boundary_cache.m](../src/+mwecmass/+geometry/precompute_boundary_cache.m),
[geometry/extract_isocurve_at_z.m](../src/+mwecmass/+geometry/extract_isocurve_at_z.m),
[geometry/extract_midplane_profile.m](../src/+mwecmass/+geometry/extract_midplane_profile.m),
[geometry/isocurve_devsurf.m](../src/+mwecmass/+geometry/isocurve_devsurf.m),
[geometry/isocurve_revsurf.m](../src/+mwecmass/+geometry/isocurve_revsurf.m), and
[geometry/isocurve_bloftsurf.m](../src/+mwecmass/+geometry/isocurve_bloftsurf.m).

## 4. Production hydrostatics: contour tables and strip integration

### 4.1 Adaptive waterplane table

The production configuration builder samples horizontal contours between the sampled hull limits. Its
initial grid count is

\[
n_c=\max\left(20,\left\lceil\frac{z_{max}-z_{min}}{\Delta z_{target}}\right\rceil\right),
\]

or 40 when no positive target spacing is supplied. It excludes the exact endpoints by `10^-4 m`,
extracts a contour at every coarse elevation, and computes:

\[
A_w=\frac12\left|\sum_i (x_i y_{i+1}-x_{i+1}y_i)\right|,
\]

\[
I_{wp,xx}=\frac1{12}\left|\sum_i q_i(y_i^2+y_i y_{i+1}+y_{i+1}^2)\right|,
\]

\[
I_{wp,yy}=\frac1{12}\left|\sum_i q_i(x_i^2+x_i x_{i+1}+x_{i+1}^2)\right|,
\quad q_i=x_i y_{i+1}-x_{i+1}y_i.
\]

Polygon perimeter `P(z)` is the sum of closed-edge lengths. Intervals whose
`|Delta A_w/Delta z|` exceeds 30% of the largest coarse-grid gradient receive three interior samples.
The builder adds exact zero-area endpoints at the hull limits.

Every strip boundary missing from this grid is subsequently evaluated as an exact contour-table node.
The augmentation appends and jointly sorts `A_w`, `I_wp,xx`, `I_wp,yy`, `V_sub`, `CB_z`, and—when
available—`S_wet`; it then rebuilds the cumulative volume and centroid tables from scratch. It does
**not** append corresponding values to the perimeter table `P(z)`. This avoids interpolating area and
moment quantities across sharp shoulders at strip edges, but leaves `P_table` on the original grid.

That length mismatch is consequential in C1 modular-precast post-processing: `extract_strip_geometry`
uses `P_table` only when its length matches the augmented `A_w` grid. For C1 it therefore estimates the
slope correction from the equivalent-area radius `sqrt(A_w/pi)` rather than from `dA_w/dz divided by P`.

### 4.2 Displaced volume and center of buoyancy

For each table elevation `z_k`, production displacement is cumulative trapezoidal integration:

\[
V_{sub}(z_k)=\int_{z_{min}}^{z_k} A_w(z)\,dz,
\qquad
CB_{z,b}(z_k)=\frac{\int_{z_{min}}^{z_k}zA_w(z)\,dz}{V_{sub}(z_k)}.
\]

The actual code uses `trapz` over all table nodes at or below `z_k`. For volume below `10^-10 m^3`, it
sets volume to zero and uses `z_min` as a finite centroid placeholder. It asserts that the rebuilt volume
table does not decrease by more than `10^-9 m^3` between nodes.

At runtime, `properties_3d` linearly interpolates `A_w`, `I_wp`, `V_sub`, `CB_z`, and wetted area at
`z_sub_top=min(-v_s,z_max)`. World-frame buoyancy height is `CB_z=CB_z,b+v_s`. This lookup is the primary
hydrostatic path used by the optimizer and realization solvers.

The table-integrated maximum displaced volume replaces the preliminary divergence-theorem volume in the
configuration. This is deliberate: for an open surface set, divergence-theorem closure can include a
virtual cap, whereas `integral A_w dz` counts the represented cross-sections.

### 4.3 Horizontal strip geometry

For strip `i=[z_i,z_{i+1}]`, `compute_strip` collects all table nodes inside the interval and inserts both
endpoints. It then evaluates

\[
V_i=\int A_w(z)dz,\qquad
\bar z_i=\frac{\int zA_w(z)dz}{V_i},\qquad
Q_{z^2,i}=\int z^2A_w(z)dz,
\]

\[
Q_{x^2,i}=\int I_{wp,yy}(z)dz,\qquad
Q_{y^2,i}=\int I_{wp,xx}(z)dz.
\]

The volume moments become

\[
I_{yy,i}^{(vol)}=Q_{x^2,i}+Q_{z^2,i},\quad
I_{xx,i}^{(vol)}=Q_{y^2,i}+Q_{z^2,i},\quad
I_{zz,i}^{(vol)}=Q_{x^2,i}+Q_{y^2,i}.
\]

These are geometric volume moments in `m^5`; density is applied later. A degenerate strip returns zero
volume and moments and the interval midpoint as its centroid.

Sources:
[driver/build_hydrostatic_tables.m](../src/+mwecmass/+driver/build_hydrostatic_tables.m),
[driver/build_strip_geometry_tables.m](../src/+mwecmass/+driver/build_strip_geometry_tables.m),
[hydrostatics/waterplane_properties.m](../src/+mwecmass/+hydrostatics/waterplane_properties.m), and
[hydrostatics/compute_strip.m](../src/+mwecmass/+hydrostatics/compute_strip.m).

## 5. Divergence theorem, Gauss–Legendre quadrature, and clipping

### 5.1 What is—and is not—the production displacement method

The surface-integral engine applies `n×n` Gauss–Legendre quadrature on `[0,1]^2`. With oriented area
vector `n dA=(S_u×S_v)du dv`, a source surface contributes

\[
V_s=\frac13\iint S\cdot(S_u\times S_v)\,du\,dv.
\]

Coordinate first moments use `x^2 n_x/2`, `y^2 n_y/2`, and `z^2 n_z/2`; second volume moments use the
corresponding cubic terms divided by three. A coarse `5×5` signed-volume calculation selects the source
orientation. Mirror contributions are generated analytically, and detected open edges can receive flat
cap contributions.

This machinery is used to obtain initial whole-hull centroid and moments, support geometry diagnostics,
compute wetted surface area, and provide a fallback when hydrostatic tables are absent. `compute_strip`
also retains a surface-integral fallback. It is **not** the primary source of optimizer-loop `V_sub`,
`CB_z`, or strip geometry when the standard configuration tables exist. Partial-surface quadrature can
poorly resolve a waterline discontinuity and can be sensitive to orientation or implicit caps; the
cross-sectional table path was implemented to avoid those failure modes.

The Gauss–Legendre helper uses the Golub–Welsch eigensystem and maps nodes from `[-1,1]` to `[0,1]`. An
`n`-point rule is exact for polynomials through degree `2n-1`; parametric hull integrands are generally
not polynomials, so quadrature order remains an accuracy control, not an exactness guarantee.

### 5.2 Horizontal and panel clipping

The code contains three related clipping implementations:

| routine | domain | operation |
|---|---|---|
| `clip_polygon_at_z` | 2-D `[x,z]` polygon | Sutherland–Hodgman-style traversal retaining `z<=cut` or `z>=cut`; inserts linear edge intersections |
| `internal/clip_z` | paired `x,z` vectors | compact equivalent used by realization/figure helpers |
| `mesh/split_panel_at_z` | 3-D triangle/quad panel | traverses four perimeter edges, retains submerged vertices, inserts and snaps crossings to `z_wl` |

Panel clipping preserves input winding. A three-vertex output is stored as `[v1 v2 v3 v3]`; polygons with
more than four vertices are fan-triangulated. The mesh routine uses a `10^-10 m` submerged-side tolerance
and a `10^-14 m` edge-height guard. The 2-D polygon routine uses `10^-9 m`, clamps intersection parameters
to `[0,1]`, deduplicates at `10^-6 m`, and returns empty when fewer than three output points remain.

For material offsets, `offset_polygon` first normalizes winding to counter-clockwise, averages adjacent
inward edge normals, and applies a bounded miter correction. When the half-angle cosine is at most 0.33,
the miter distance is capped at three times the requested offset. Offset validity is checked afterward;
a collapsed, inverted, or larger-than-outer inner polygon is treated as locally solid.

Sources:
[hydrostatics/surface_integral.m](../src/+mwecmass/+hydrostatics/surface_integral.m),
[hydrostatics/compute_hull.m](../src/+mwecmass/+hydrostatics/compute_hull.m),
[internal/gauss_legendre.m](../src/+mwecmass/+internal/gauss_legendre.m),
[geometry/clip_polygon_at_z.m](../src/+mwecmass/+geometry/clip_polygon_at_z.m),
[mesh/split_panel_at_z.m](../src/+mwecmass/+mesh/split_panel_at_z.m), and
[internal/offset_polygon.m](../src/+mwecmass/+internal/offset_polygon.m).

## 6. Mass properties and flotation

### 6.1 Strip mass, center of gravity, and inertia

The optimized density is piecewise constant by horizontal strip. Densities are clipped to their configured
bounds inside the property evaluator as a defensive guard; normally the optimizer bounds already enforce
the same interval. For precomputed geometric properties,

\[
m_i=\rho_i V_i,\qquad M=\sum_i m_i,
\qquad CG_{z,b}=\frac{\sum_i m_i\bar z_i}{M}.
\]

The exact algebraic mass closure of this model is therefore

\[
M-\sum_i\rho_iV_i=0
\]

up to floating-point summation. There is no independent “unassigned” mass term. World-frame
`CG_z=CG_z,b+v_s`; symmetry assumptions set `CG_x=CG_y=0` in the 3-D evaluator.

Density-weighted origin moments are

\[
I_{yy,O}=\sum_i\rho_i I_{yy,i}^{(vol)},\quad
I_{xx,O}=\sum_i\rho_i I_{xx,i}^{(vol)},\quad
I_{zz,O}=\sum_i\rho_i I_{zz,i}^{(vol)}.
\]

Because the modeled CG lies on the `z` axis,

\[
I_{xx,G}=\max(0,I_{xx,O}-M CG_{z,b}^2),\quad
I_{yy,G}=\max(0,I_{yy,O}-M CG_{z,b}^2),\quad
I_{zz,G}=\max(0,I_{zz,O}).
\]

The `max(0,...)` operations clamp every negative result, regardless of magnitude; they are not gated by a
round-off threshold. If configured CAD inertia values differ beyond the permitted percentage, the
evaluator warns and may replace the computed diagonal inertias with those overrides.

### 6.2 Free-floating equilibrium

The displaced-water mass is

\[
M_b=\rho_w V_{sub}(-v_s).
\]

Stage 2 and both material realization solvers enforce the dimensionless equality

\[
c_{eq}=\frac{M}{M_b}-1=0.
\]

This is a nonlinear solver constraint. “Exact mass closure” refers to the model identities above; a
reported optimized point satisfies flotation only to the solver's finite constraint tolerance. The code
also stores the signed discrepancy `M-M_b` and, in realization, its absolute value.

There is no mooring or PTO stiffness in the shipped free-floating model:

\[
K_{PTO}=0,\qquad K_{total}=K_{hydro}.
\]

Consequently surge is unrestored and its natural period is infinite.

Source: [hydrostatics/properties_3d.m](../src/+mwecmass/+hydrostatics/properties_3d.m).

## 7. Stability, stiffness, hydrodynamics, and periods

### 7.1 `CB`, `BM`, `KM`, and `GM`

For longitudinal stability in pitch,

\[
BM_L=\frac{I_{wp,yy}}{V_{sub}},\qquad
KM=CB_z+BM_L,\qquad
GM_L=KM-CG_z.
\]

All elevations in these expressions are in the same world frame. Translation by `v_s` shifts `CB_z` and
`CG_z` equally and therefore cancels directly from their difference, while the selected submerged geometry
and its waterplane still change with `v_s`.

> **Axis warning:** `I_wp,yy` is calculated about body `x=0`, not shifted to the instantaneous center of
> flotation. The expression equals a centroidal waterplane moment only when the waterplane centroid lies on
> `x=0`, as it does for the axisymmetric C1 hull.

### 7.2 Hydrostatic stiffness

For the `[surge,heave,pitch]` subset,

\[
K_{hydro}=\operatorname{diag}(0,K_{33},K_{55}),
\]

\[
K_{33}=\rho_w g A_w,
\qquad
K_{55}=\begin{cases}MgGM_L,&GM_L>0,\\0,&GM_L\le0.\end{cases}
\]

At flotation equilibrium, `M=rho_w V_sub`; away from equilibrium during optimization, the implementation
deliberately uses physical candidate mass `M` in `K55`.

### 7.3 Hydrodynamic cache interpolation and CG transfer

The cache supplies infinite-frequency added mass `A∞` and band-averaged radiation damping `Bbar` at a set
of `vertical_shift` nodes. The code clamps the query to the cached interval, linearly interpolates every
matrix element, and symmetrizes full interpolated matrices. Before the optional CG transfer, it clamps the
three scalar added-mass values and the diagonals of both full matrices nonnegative. A one-node cache is
used directly. Missing full matrices fall back to diagonal `A∞` and zero `Bbar`.

Cached coefficients are referenced to the BEM/HAMS center of gravity. For target vertical CG offset
`Delta z`, generalized velocities obey

\[
q_O=Tq_G,\qquad
T=\begin{bmatrix}1&0&-\Delta z\\0&1&0\\0&0&1\end{bmatrix},
\]

so energy invariance gives the congruence transfer

\[
A_G=T^T A_O T,\qquad B_G=T^T B_O T.
\]

The transfer routine symmetrizes matrices before and after the congruence. After a nontrivial transfer,
`interpolate_at_draft` re-extracts `A11`, `A33`, and `A55` from the transformed matrix and clamps those
three scalar outputs nonnegative; it does not repeat the full-matrix diagonal clamp after the transform.
If any operation in `interpolate_at_draft` escapes its internal helpers and reaches its outer catch, the
function warns and returns zero scalar coefficients plus zero `A` and `B` matrices. `A∞` enters the
natural-period mass matrix; `Bbar` is exported and plotted but does not enter the undamped eigenproblem
described below.

Sources:
[bem/interpolate_at_draft.m](../src/+mwecmass/+bem/interpolate_at_draft.m),
[bem/interpolate_matrix.m](../src/+mwecmass/+bem/interpolate_matrix.m), and
[bem/transform_to_cg.m](../src/+mwecmass/+bem/transform_to_cg.m).

### 7.4 Uncoupled and coupled natural periods

The logged per-axis comparisons are

\[
T_h^{unc}=2\pi\sqrt{\frac{M+A_{33}^{\infty}}{K_{33}}},\qquad
T_p^{unc}=2\pi\sqrt{\frac{I_{yy,G}+A_{55}^{\infty}}{K_{55}}}.
\]

A stiffness at or below `10^-6` produces `Inf`. These values are retained as
`heave_uncoupled` and `pitch_uncoupled`.

The periods consumed by the Stage 2 objective are the coupled undamped modes. The engine solves

\[
K\phi_j=\lambda_j M_v\phi_j,
\qquad
M_v=\operatorname{diag}(M,M,I_{yy,G})+A^\infty,
\qquad
T_j=\frac{2\pi}{\sqrt{\lambda_j}}.
\]

Nonpositive eigenvalues map to `Inf`. Modes are labeled using each coordinate's signed kinetic-energy
share,

\[
s_{ij}=100\frac{\phi_{ij}(M_v\phi_j)_i}{\phi_j^T M_v\phi_j}.
\]

Every column sums to 100%; individual entries can be negative because an off-diagonal mass term is split
between its coupled coordinates. The largest share names the mode. If a coordinate with positive diagonal
stiffness has no mode with at least 50% share, the code raises `AmbiguousModeLabel` rather than silently
assigning it. Since surge has zero stiffness, its mode remains unrestored.

Source:
[hydrostatics/coupled_periods_by_share.m](../src/+mwecmass/+hydrostatics/coupled_periods_by_share.m).

## 8. Optimization formulation

### 8.1 Stage 1 modes

`stage1_mode` selects one of four implemented warm-start strategies:

| mode | executed method | status |
|---|---|---|
| `sweep` | full 3-D property screen at cached draft nodes, followed by short fixed-draft SQP on the best candidates | shipped default and C1 reference path |
| `skip` | pass the configured initial shift and densities directly to Stage 2 | cold-start option |
| `oneshot` | one 2-D surrogate SQP followed by a 3-D comparison | diagnostic/alternative warm start |
| `trained` | repeated 2-D surrogate SQP with PID updates to volume and CG correction factors, each iteration checked in 3-D | alternative calibrated-surrogate path |

The `sweep` banner retains historical “2-D density search” wording, but its decisions are made with
`properties_3d`. The 2-D evaluator called at the end only populates reporting fields.

For each sweep node, Tier 1 selects one density that attempts flotation balance:

\[
\rho_{uniform}=\frac{\rho_wV_{sub}}{\sum_i V_i}.
\]

In constructability mode, the wall strip is fixed at hull-material density and the remaining required mass
is spread uniformly across the other strips. Bounds can clip these analytical densities. Tier 1 marks a
candidate feasible when relative mass error is below 5% and `GM_L>=gm_min`. Candidates are ranked by
`|T_h-T_h*|`; the configured best count receives a fixed-`v_s` SQP refinement using the Stage 2 objective
and constraints. A refined point is accepted only when it improves the objective; its feasibility flag
also requires positive exit status, relative mass error below 1%, and the GM floor.

`trained` first estimates multiplicative corrections `k_vol` from the 3-D/2-D displaced-volume ratio and
`k_gm` from the 3-D/2-D CG ratio. It then alternates 2-D SQP, optional hydrodynamic enrichment, 3-D
validation, and bounded PID updates. It converges only with a positive 2-D exit status and either both
error tolerances satisfied or stable correction factors after the warm-up interval. Three consecutive
iterations with both factors saturated cause an early best-effort exit. With live BEM disabled, validation
interpolates the existing cache; it does not generate exact new coefficients.

The 2-D model is a midplane-extrusion surrogate with an area-matched effective transverse width, finite
horizontal strips, `k_vol` volume scaling, and a `k_gm` CG bias. It is not the production 3-D hydrostatic
model. Its first design variable is `vertical_shift`, not draft.

Sources:
[optim/stage1_sweep.m](../src/+mwecmass/+optim/stage1_sweep.m),
[optim/stage1_screen_draft.m](../src/+mwecmass/+optim/stage1_screen_draft.m),
[optim/stage1_oneshot.m](../src/+mwecmass/+optim/stage1_oneshot.m),
[optim/stage1_trained.m](../src/+mwecmass/+optim/stage1_trained.m), and
[hydrostatics/properties_2d.m](../src/+mwecmass/+hydrostatics/properties_2d.m).

### 8.2 Stage 2 variables and bounds

Stage 2 optimizes `x=[v_s,rho_1,...,rho_N]`. When automatic bounds are selected, the builder sets
`v_s,min=-z_max+0.1 m` and `v_s,max=-z_min-0.1 m`, keeping the body-frame waterline 0.1 m inside the
hull limits. A nonempty manual `vertical_shift_bounds` input is accepted verbatim: the builder does not
reapply that margin or validate physical immersion, so the caller owns its suitability. Every density
stays inside the configured ballast-density interval. Constructability can raise individual lower bounds;
its structural wall strip is fixed by making its lower and upper bounds equal to the hull-material density.

For a constructable non-wall strip, the lower bound approximates a minimum-thickness normal-offset shell.
Using the minimum radius from the cross-section centroid to its boundary as an equivalent profile
`r(z)`, the builder uses

\[
r_{in}(z)=\max\left(0,r(z)-t_{min}\sqrt{1+r'(z)^2}\right),
\]

\[
V_{shell}=\pi\int\left(r^2-r_{in}^2\right)dz,
\qquad
\rho_{min}=\frac{\rho_{hull}V_{shell}+\rho_{void}(V_i-V_{shell})}{V_i}.
\]

`rho_min` is clamped between constituent densities and combined with the global ballast lower bound.
Where the normal offset reaches the axis, the slice is solid. For a non-axisymmetric cross-section,
the centroid-to-boundary minimum-radius substitution is an approximation; the implementation describes
the resulting shell volume as a lower-bound construction model. The minimum constructable mass used by
the hard constraint is `sum(V_i rho_min,i)` with the wall strip pinned.

If the Stage 1 warm start violates the GM floor, a Phase A SQP first holds `v_s` fixed and adjusts only the
density distribution. Phase A is accepted only on positive exit and a satisfied GM floor. The main solve
then uses the configured `fmincon` algorithm (SQP by default), scaled variables, a `10^-8` constraint and
optimality tolerance, and a `10^-10` step tolerance.

### 8.3 Objective and range penalty

For `q in {GM_L,T_h,T_p}`, define

\[
r_q=\frac{q-q^*}{\tfrac12(q_{hi}-q_{lo})}.
\]

The dimensionless objective is

\[
J=\phi(r_{GM})+\phi(r_h)+\phi(r_p),
\]

with the `C1`-continuous penalty

\[
\phi(r)=
\begin{cases}
r^2,&|r|\le1,\\
1+2\delta+k_{amp}\delta^2,&|r|>1,\quad\delta=|r|-1.
\end{cases}
\]

Value and slope match at `|r|=1`; `k_amp` increases only the curvature beyond the target band. The Stage 2
guard returns the configured penalty when `GM_L` is `NaN`, or when either heave or pitch period is `NaN`
or infinite. It does not explicitly guard an infinite `GM_L`; that value propagates through the penalty
arithmetic instead.

The GM and period ranges are **soft objective bands**, not hard constraints. The only period acceptance
check is a reported quality metric after optimization; it does not determine the Stage 2 `converged`
boolean.

### 8.4 Hard constraints

`fmincon` receives inequalities `c<=0` and equality `ceq=0`:

\[
c_{GM}=1-\frac{GM_L}{GM_{min}},
\]

\[
c_{ratio,i}=\frac{\rho_i}{\rho_{i+1}+1}-R_{max},\qquad
c_{mono,i}=\frac{\rho_{i+1}-\rho_i}{\rho_{max}},
\]

\[
c_{eq}=\frac{M}{\rho_wV_{sub}}-1.
\]

Adjacent constraints exclude any pair touching a fixed constructability wall strip. When a positive
minimum constructable mass exists, the additional inequality is

\[
c_{minmass}=1-\frac{M}{M_{min,constructable}}.
\]

The `+1 kg/m^3` in the adjacent ratio denominator is an implemented smoothness guard. Density monotonicity
requires non-increasing density upward.

The constraint function has an outer catch that returns fixed positive violations for exceptions raised
through it. That is not a guarantee for every invalid candidate: `properties_3d` catches many internal
errors itself and returns defaults, after which `M/M_b-1` can evaluate to `NaN` (for example, `0/0`)
without triggering the outer catch. Solver logs and saved constraint residuals remain authoritative.

### 8.5 Convergence semantics

The main solver exit flag alone is not the saved Stage 2 convergence decision. After the solve,
`check_3d_convergence` requires all of:

- exact-sign monotonicity at the reported point (`diff(rho)<=0`, excluding wall-adjacent pairs);
- absolute flotation discrepancy below 10 kg;
- `GM_L>gm_min` (strictly greater in this post-check);
- either exit flag 1, or exit flag 0/2 with constraint violation below `10^-6` and first-order optimality
  below `10^-2`.

Period errors below 10% are stored as `periods_acceptable`, but are not part of the logical conjunction.
Therefore a run can be marked converged while a soft period target is missed, or marked not converged at a
numerically useful boundary point. Consumers should inspect the saved exit flag, quality metrics, residuals,
and objective together.

Sources:
[optim/stage2_bounds.m](../src/+mwecmass/+optim/stage2_bounds.m),
[optim/stage2_objective.m](../src/+mwecmass/+optim/stage2_objective.m),
[optim/range_penalty.m](../src/+mwecmass/+optim/range_penalty.m),
[optim/stage2_constraints.m](../src/+mwecmass/+optim/stage2_constraints.m), and
[optim/check_3d_convergence.m](../src/+mwecmass/+optim/check_3d_convergence.m).

## 9. Realization models

### 9.1 Preliminary

`preliminary` performs no material solve. It returns the Stage 2 optimized equivalent-density profile and
its `Final3D` properties unchanged. Density values in this mode are mathematical bulk densities, not a
fabrication prescription.

Source: [realise/preliminary/run.m](../src/+mwecmass/+realise/+preliminary/run.m).

### 9.2 Thin shell with ballast level

The thin-shell solve uses three design variables:

\[
x_r=[v_s,t,z_{ballast}]^T.
\]

It minimizes only heave and pitch range penalties, subject to flotation equality and the GM floor. GM is
not an objective term in this realization stage. Thickness is bounded below by the fabrication minimum
and above by 95% of a radius-derived limit. `z_ballast` lies just inside the sampled hull limits.

The inner realization SQP evaluates the **uncoupled** per-axis formulas for its heave and pitch objective.
After realization, the shared property builder recomputes and exports coupled, kinetic-share-labeled
periods at the realized mass, CG, and draft. The inner objective periods and final exported periods can
therefore differ when surge-pitch added-mass coupling is nonzero.

The SQP flotation equality and GM inequality define the intended formulation, but they are not the gate
used to decide whether the returned thin-shell point replaces `final_props`. The saved `feasible` flag is
only `isfinite(GM) && isfinite(T_heave) && isfinite(T_pitch)`. It does not require a positive `exitflag`,
small flotation residual, or satisfaction of the GM floor. Consequently, a finite point from a failed or
nonconverged SQP can replace the optimizer properties. Users must inspect `steel_data.exitflag`,
`steel_data.mass_balance_error_pct`, `steel_data.GM_realised`, and `steel_data.residuals` before treating
the realized result as accepted.

#### Inner boundary and slope correction

At each geometry-grid elevation, the outer contour is offset inward. The implementation estimates

\[
\frac{dr}{dz}\approx\frac{|dA_w/dz|}{P},\qquad
\cos\alpha=\frac1{\sqrt{1+(dr/dz)^2}},\qquad
t_{in\mbox{-}plane}=\frac{t}{\cos\alpha}.
\]

`dA_w/dz` uses a `0.01 m` finite difference. `cos(alpha)` is floored by the configured maximum slope
factor; a near-top cap uses 0.1 before that bound. The inward contour uses the miter offset of Section 5.
If the tabulated outer area is positive but contour extraction returns fewer than three points,
`inner_properties_at_z` raises `NoIsocurveAtZ`; the solver does not silently convert a missing outer
contour into material. For a successfully extracted outer contour, an offset with fewer than three points,
inner area at or below `10^-10 m^2`, or inner area greater than or equal to the outer area becomes a solid
slice (`A_inner=I_inner=0`). A nonfinite inner moment, by contrast, falls back to the corresponding outer
moment multiplied by `(A_inner/A_outer)^2`. Inner moments are then capped at the outer moments.

This is a contour-offset approximation. The relation `dA/dz≈P dr/dz` is exact only for restricted shape
families and uniform normal motion; corners and rapidly changing non-axisymmetric sections are approximate.

#### `z_ballast` partition and two solid densities

`z_ballast` is inserted as an exact interpolation knot in the geometry grid. Below it, the full outer section
is solid ballast. Above it, the material occupies the annulus and the inner section is air:

\[
V_{ballast}=\int_{z_{min}}^{z_{ballast}}A_o dz,
\]

\[
V_{shell}=\int_{z_{ballast}}^{z_{max}}(A_o-A_i)dz,
\qquad
V_{air}=\int_{z_{ballast}}^{z_{max}}A_i dz.
\]

Mass closes by region:

\[
M=\rho_{ballast}V_{ballast}+\rho_{shell}V_{shell}+\rho_{air}V_{air}.
\]

The source preserves an equivalent aggregate expression when `rho_ballast=rho_shell` so the equal-density
default follows the original arithmetic order. First moments and all three diagonal inertias are
superposed from the same region integrals, then transferred from the body origin to the CG. The solver
requires `rho_air<min(rho_shell,rho_ballast)` and `rho_ballast>=rho_shell`; the latter guarantees monotonicity of
the generalized analytical ballast-seed integral. For unequal solid densities, the seed inverts the exact
quadratic cumulative on each interval of a piecewise-linear grid integrand. End clamps are bounded
fallbacks when the requested mass lies outside the seed's achievable range.

The final realized properties are rebuilt at the realized `v_s` and CG, including a fresh hydrodynamic
interpolation and coupled-period calculation.

Sources:
[realise/thin_shell/inner_properties_at_z.m](../src/+mwecmass/+realise/+thin_shell/inner_properties_at_z.m),
[realise/thin_shell/hull_slope_cos_at_z.m](../src/+mwecmass/+realise/+thin_shell/hull_slope_cos_at_z.m),
[realise/thin_shell/integrate_split.m](../src/+mwecmass/+realise/+thin_shell/integrate_split.m),
[realise/thin_shell/evaluate_design_point.m](../src/+mwecmass/+realise/+thin_shell/evaluate_design_point.m), and
[realise/thin_shell/solve.m](../src/+mwecmass/+realise/+thin_shell/solve.m).

### 9.3 Modular precast UHPC

The modular-precast realization solves

\[
x_r=[v_s,z_{ballast},t_1,\ldots,t_{N_{nonwall}}]^T,
\]

where each non-wall strip has an offset thickness. The designated structural wall strip is permanently
solid and is excluded from thickness variables. As in thin shell, the objective contains heave and pitch
penalties; flotation and the GM floor are hard constraints.

As in the thin-shell branch, the modular-precast SQP objective uses uncoupled per-axis periods. Final
reported properties are rebuilt through the shared coupled-mode calculation.

The same limited `feasible` test is used here: finite GM, heave period, and pitch period. It does not test
the modular SQP `exitflag`, flotation residual, or GM constraint. `modular_precast.run` passes a finite
result through the shared realized-property builder even after a failed/nonconverged solve; only a false
`feasible` flag causes that builder to retain the optimizer properties. Inspect
`constructability.exitflag`, `constructability.mass_balance_error_pct`,
`constructability.GM_realised`, and `constructability.residuals` before accepting final realized values.

At each elevation sample, the material classification is:

| condition | section physics |
|---|---|
| designated wall strip | fully solid UHPC |
| strip explicitly promoted solid | fully solid UHPC |
| `z<=z_ballast` | fully solid UHPC |
| otherwise, valid inward offset | UHPC annulus plus air/void interior |
| solver grid: collapsed, tiny, or too-large inner polygon | fully solid fallback |

During the constrained solve, a positive tabulated outer area with a missing extracted outer contour raises
the same `NoIsocurveAtZ` error as thin shell. The more permissive missing-contour behavior occurs only in
the later `extract_strip_geometry` pass: an unavailable or invalid sampled contour is assigned no void and
falls back to the available outer-section material representation for post-solve strip reporting. That
post-processing fallback must not be mistaken for the solver's behavior.

The solver's grid uses the same slope-corrected inner-contour machinery as thin shell. After optimization,
non-wall strips lying completely below `z_ballast` are promoted to solid for downstream extraction.

For each extracted strip and sample,

\[
A_{UHPC}=A_o-A_i,\qquad A_{void}=A_i
\]

for an annular section, or `A_UHPC=A_o`, `A_void=0` for a solid one. Trapezoidal integration gives
`V_UHPC` and `V_void`, and

\[
M_i=\rho_{UHPC}V_{UHPC,i}+\rho_{void}V_{void,i},
\]

\[
\rho_{eff,i}=\frac{M_i}{V_{UHPC,i}+V_{void,i}}.
\]

Using the same denominator integrator as both component volumes guarantees that `rho_eff` remains a
mass-weighted mixture between the two constituent densities. The reported void contour is the same
slope-corrected contour used in the mass calculation.

> **Authority warning:** the constrained global solve is the source of final mass, GM, draft, and periods.
> The later per-strip extraction uses a separate sampled integration for visualization and diagnostics.
> Its sum is compared with the authoritative solve and may differ; the code labels a difference below 2%
> “OK” and explicitly states that strip data is visual only. Do not replace final properties with the
> diagnostic strip sum.

The pre-solve constructability guard estimates a minimum mass using the wall-solid/minimum-jacket model
and a maximum mass from fully solid UHPC. Targets beyond tolerance around this achievable bracket raise
an error before SQP.

Sources:
[realise/modular_precast/solve_and_extract.m](../src/+mwecmass/+realise/+modular_precast/solve_and_extract.m),
[realise/modular_precast/build_geometry_grid.m](../src/+mwecmass/+realise/+modular_precast/build_geometry_grid.m),
[realise/modular_precast/evaluate_design_point.m](../src/+mwecmass/+realise/+modular_precast/evaluate_design_point.m),
[realise/modular_precast/solve.m](../src/+mwecmass/+realise/+modular_precast/solve.m), and
[realise/modular_precast/extract_strip_geometry.m](../src/+mwecmass/+realise/+modular_precast/extract_strip_geometry.m).

## 10. Numerical safeguards and failure behavior

The principal guards are part of the implemented method, not incidental logging:

| location | guard / fallback |
|---|---|
| waterplane polygon | fewer than three points returns zero area and moments |
| cumulative hydrostatics | tiny volume returns zero and a finite centroid placeholder; monotonicity checked |
| strip table | exact endpoints inserted; negative interpolated areas/moments clamped to zero |
| inertia | all negative CG moments clamped to zero, without a magnitude threshold |
| hydro cache | query clamped; matrices symmetrized; pre-transfer diagonals clamped; outer failure returns zero matrices |
| dynamics | near-zero stiffness gives infinite period; ambiguous restored-mode label raises an error |
| optimizer objective | `NaN` GM or nonfinite heave/pitch period returns a large guard penalty |
| optimizer constraints | thrown exceptions return fixed violations; internally caught defaults can still yield `NaN` equality residuals |
| material offset | missing outer contour errors in solvers; collapsed/tiny/too-large inner polygon becomes solid |
| realization grid | more than 5% degenerate sections raises an error |

Several geometry helpers catch an exception and return an empty or zero result. That behavior prevents a
single malformed polygon from crashing low-level plotting or fallback code, but in production it may
surface later as zero waterplane, an objective guard penalty, or infeasibility. Review warnings and the
saved residuals; a numerically completed run is not automatically a physically valid one.

## 11. Assumptions and limitations

1. **SI geometry is required.** The parser records `Units:` but does not scale coordinates.
2. **Symmetric-body mass model.** The 3-D optimizer forces `CG_x=CG_y=0` and retains only diagonal rigid-
   body inertia. It is not a general asymmetric mass-distribution solver.
3. **Origin-centered waterplane inertia.** `BM_L` is exact as implemented for waterplanes centered on
   `x=0`; an off-center waterplane needs a center-of-flotation parallel-axis correction that is not present.
4. **Contour topology.** Angular sorting assumes a single well-behaved outer loop. Holes, disconnected
   loops, and strongly concave sections can be misrepresented.
5. **First-crossing iso-curves.** Non-monotone revolution and axial-loft profiles can have multiple roots;
   the implemented helpers select the first.
6. **Sampled geometry.** Surface extents, profiles, contour tables, and realization grids have finite
   resolution. Adaptive `z` refinement responds only to waterplane-area gradient.
7. **MultiSurf subset.** Only the entity and surface types listed in Section 3 are evaluatable. Unknown
   constructs are not preserved as generic geometry.
8. **Developed and lofted surfaces are approximations.** `DevSurf` is a ruled same-parameter blend;
   `BLoftSurf` synthesizes its knot vector and assumes corresponding section parameters.
9. **Free-floating, undamped natural periods.** PTO/mooring stiffness is zero, and radiation damping is not
   included in the eigenproblem. The results are undamped hydrostatic/added-mass periods.
10. **Infinite-frequency added mass.** Period calculations use `A∞`, not frequency-dependent added mass
    evaluated at the resulting natural frequency.
11. **Soft target ranges.** GM target and period bands shape the objective; only the GM floor and the
    constraints in Section 8.4 are hard requirements.
12. **Local nonlinear optimization.** SQP results depend on bounds, scaling, tolerances, and warm start.
    Density design vectors should not be assumed unique in flat objective directions.
13. **Realization approximations.** Mitered polygon offsets and the `dA/dz` slope correction approximate a
    constant normal shell thickness. Collapsed offsets are intentionally converted to solid material.
14. **Diagnostic strip reconstruction.** Modular-precast post-solve strip sums are not the authoritative
    final mass properties.

## 12. Source map for the implemented method

| concern | primary implementation |
|---|---|
| MultiSurf parsing and surface evaluation | [MS2Parser.m](../src/+mwecmass/+geometry/MS2Parser.m) |
| production contour/table hydrostatics | [build_hydrostatic_tables.m](../src/+mwecmass/+driver/build_hydrostatic_tables.m) |
| strip geometry | [build_strip_geometry_tables.m](../src/+mwecmass/+driver/build_strip_geometry_tables.m), [compute_strip.m](../src/+mwecmass/+hydrostatics/compute_strip.m) |
| candidate mass, stability, and periods | [properties_3d.m](../src/+mwecmass/+hydrostatics/properties_3d.m) |
| hydrodynamic interpolation and transfer | [interpolate_at_draft.m](../src/+mwecmass/+bem/interpolate_at_draft.m), [transform_to_cg.m](../src/+mwecmass/+bem/transform_to_cg.m) |
| Stage 1 and Stage 2 | [optim/run.m](../src/+mwecmass/+optim/run.m) |
| shared objective penalty | [range_penalty.m](../src/+mwecmass/+optim/range_penalty.m) |
| realized property rebuild | [build_realised_properties.m](../src/+mwecmass/+realise/build_realised_properties.m) |
| constructability density floor | [compute_perpendicular_shell_volume.m](../src/+mwecmass/+hydrostatics/compute_perpendicular_shell_volume.m) |

This source map intentionally names the executable modules that define the method. It does not imply that
every helper has been independently validated for every geometry family.
