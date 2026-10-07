# Interface contract for the exact-geometry, Stage-3 and STEP tasks

> Binding for every task of `2026-10-06-exact-geometry-stage3-step.md` from the lanes on; read
> `AGENTS.md` first. Built on the merged T1 (`+solid/outer_rows.m`, `surface_normals.m`; used by
> T2a, T2b, T3 tests) and T9 (`+step/write_step.m`, `validate_brep.m`, `tests/step_check.py`; used by T3
> tests, T10). The contract wins where a plan task section differs (§3). Cited S1–S8, F1–F14, I1–I9.

## 0. Conventions

- **Units:** SI, lengths in metres, everywhere (structs, tables, STEP). Densities kg/m³.
- **Body frame:** the `.ms2` coordinates (C1: keel z = −3.25, top z = 1.10). **World frame:** the
  same x, y; z_world = z_body + vs, vs = `vertical_shift` = x(1); z_world = 0 is the still-water
  line, so the waterline is z_body = −vs and draft = −(z_min + vs). Every field name or comment
  states its frame; geometry structs are body frame, `CG_total`, `CB`, `KM` are world frame.
- **Names (rule 10):** regions are `uhpc`, `air` (modular precast) and `ballast`, `shell`, `air`
  (thin shell); densities `rho_uhpc`, `rho_air`, `rho_ballast`, `rho_shell`; ballast top
  `z_ballast`; never `z_fill`, `rho_fill`, `rho_void`, `rho_steel` or `V_steel` for UHPC.
- **Offset distance:** a shell of design thickness t is built at d = t + eps_fit/2, eps_fit =
  0.01·t_min (AGENTS §5 item 9.2); every closed form, oracle and stand-in below uses d.
- **Errors:** `mwecmass:solid|realise|step:<Name>`. **Tests:** `tests/<area>/test_<name>.m`;
  exact-by-construction assertions only, approximation errors printed (rule 5).

## 1. Data structures

**S1 outer patch** (`geo.outer(k)`, one per entry of `model.visible_surfs`, same order, except that F1
replaces a patch in place by its pieces, consecutive, each with the patch's `name` and `source`: an
exact patch with a constant-z interval strictly inside its z range by its pieces at that interval's
ends (F1); a patch on the general path (§8) by its fitted faces; and, once any patch takes the general
path, every exact patch by its pieces at the band ends (F3b, exact; §8). Consumers map entries to
visible surfaces through `visible`, never through k. C1 and both SK fixtures: one entry per visible
surface (C1: 8)):
`name`, `source` (ultimate source entity), `type`, `flips` (`cache.mirrors(m).effective_flips`);
`surf`: T9's `'bspline'` surface as is (`degree [du dv]`, `ctrl [nu x nv x 3]` (m, body), `knots
{ku, kv}` clamped and **not renormalised** after a split (a piece keeps its parent's parameter),
`weights` `[]` or `[nu x nv]`); `outward`: S_u × S_v points out of the hull (mirrors reverse it);
`exact`: `surf` is the `.ms2` entity converted exactly; `z_of_u`: every control row has one z
(bitwise), weights separable (§5 I3); `u_range [u0 u1]`; `z_range [z(u0) z(u1)]` (ordered by u,
not by z); `offset_kind` (`rev_z` | `ruled_parallel` | `''`; F1, used by F2); `pole [1x2]`: the row
at u0 / u1 has all control points bitwise equal; `c0_u`, `c0_v`: interior knot values of full
multiplicity (du, dv), set by F1 (F2 for S2); `seam_u0`, `seam_u1`, `seam_v0`, `seam_v1`: `[patch
boundary]` of the neighbour on that boundary, empty for a pole, `[0 j]` for a row on flat region j (S1b
`flat`, §8; only on a hull with a flat region); **boundaries are numbered as the
parser's EdgeSnake** (`MS2Parser.eval_edge_snake`): 1 = v0, 2 = u1, 3 = v1, 4 = u0, of the patch as
stored (for a patch with `swap_uv`, the parser's edge e is S1 boundary 5 − e). Bitwise
equality: element-wise IEEE `==` (`isequal`), so −0 = +0. C1 values: §3, C1 oracles. Added fields
(an absent field means its default): `visible`, the index into `model.visible_surfs` of the patch
the entry comes from (default k, which holds wherever there is one entry per visible surface); in an
S2 set, the `visible` of the outer entry the piece is offset from, and for a crease fan face or a
vertex face (F2) the smallest `visible` of the outer entries that meet at its crease or vertex
(default `[]`); `fit`, for a fitted outer face (`exact` false in `geo.outer`) the S2r fields of that
face plus `dev_max` [m], the largest distance of the face from the exact surface at its check points
along the exact normal (default `[]`; always `[]` in an S2 set, whose fit report is S2r);
`swap_uv` (default false), true when F1 exchanged the entity's u and v so that z depends on u (F1):
the parser's (u, v) of a point is then the patch's (v, u).

**S1b hull geometry** `geo`: `hull_name` (deck stem), `outer` (S1 array), `z_range [z_min z_max]`
(body), `analytic` (`[]` for a real deck; stand-ins only, §3). Added field (default empty): `flat`,
the flat regions of the general path (§8: one per connected constant-z area at one height, merged
across seams and mirror planes), each `z` [m, body], `normal_z` (+1 or −1, the outward normal, out
of the hull) and `visible` (row vector: the `model.visible_surfs` indices of every patch it covers,
mirrors included); only rows of lateral faces that name it (`[0 j]`) bound it. Once per deck (and
per t_min when a patch is fitted, §8); `config.hull_solid`
(T4b, inside T0d's cached block).

**S2 inner set** (`inner`, one per distinct t): `t`, `d` = t + eps_fit/2, `eps_fit` [m], `z_range`
(body: heights of the **inner** surface the set covers, the full hollow range, §7), `z_lo` (lowest
point of the inner surface), `refit` (true when refitted on fixed knots, F2 `opts.knots_from`),
`patches` (S1 array; `outward` true when the normal points **into the void**, out of the solid;
`exact` false; crease pieces are separate patches, `seam_*` filled), `flat` (added, default empty:
constant-z regions of the inner surface as S1b `flat`, `normal_z` into the void), `report` (S2r). **S2r fit
report:** per inner patch `n_nodes`, `n_knots [u v]`, `n_passes`, `n_removed`, `n_check`,
`t_local_min`, `t_local_max` [m], `M1`, `M2`, `M3`, `M3_reason`; totals `ok`, `cap_reached`. Check
points: the F6 quadrature nodes and the midpoint of every knot span (where the integrals sample).
t_local is measured between the faces as written: outer as written (exact or fitted), inner as fitted (§8).

**S3 design** (input of F5): `mode` (`modular_precast` | `thin_shell`), `edges [N+1 x 1]` module
edges (body; `config.strip_edges`), `vs` [m], `t [N x 1]` [m] (thin shell: all equal; NaN for a
module with no void), `z_ballast` [m, body] (≤ z_min: no ballast; one value, may cross module
edges), `solid_modules` (indices solid regardless of `z_ballast`: the UHPC wall module; `[]` thin
shell). Module i layout (rule 11): full outer section solid below `z_ballast`; above it shell t(i)
+ air; `solid_modules` solid throughout; no slab at module joints. If `z_ballast` ≤ the inner
`z_lo`, the inner set is not cut: the solid (precast uhpc; thin shell ballast below `z_ballast`,
shell above) fills the section up to the inner patches, which alone close the void (C1: keel pole).

**S4 body** (output of F5)
- `design` (S3), `inner_t` (t of each S2 set used), `planes` (edges and `z_ballast` in (z_min, z_max)).
- `brep`: exactly the T9 writer struct (`vertices`, `curves`, `edges`, `surfaces`, `faces`,
  `bodies`; `uncertainty` unset = 1e-7 m). Extra face fields, which `write_step` and `validate_brep`
  ignore: `role` (`outer` | `inner` | `cap` | `ballast_top` | `joint_step`), `module`, `inside`,
  `outside` (region against / along the face normal = surface normal × (−1 if not `same_sense`);
  `exterior` outside the hull). Lateral faces (outer and inner): S1/S2 pieces cut at `planes` by
  knot insertion and at their `c0_u`, `c0_v` rows (an existing knot: no control point moves; a
  `c0_u` row of a `z_of_u` patch has one z, so its edge is shared under I2). Reason: OpenCASCADE
  `getMass` is wrong on a face with an interior full-multiplicity knot row (T9 `step_check.py`: up
  to +2.2 %; split faces −2.5e-10). Caps: `plane` faces inside (z_min, z_max), annuli with a
  hole loop. A collapsed row is omitted from the loop (three edges, or two when both rows collapse).
  The hull's own bottom and top are outer patches (C1: poles; box: horizontal faces), never caps.
- Modular precast faces: `cap` at every module joint (the full section where both sides are solid,
  else the annulus between the outer section and the inner section of smaller t, plus a uhpc/air
  disk where only one side is hollow, e.g. the void top under the solid wall module); `joint_step`
  (annulus between the two inner sections where t differs, uhpc/air); `ballast_top` is only the
  void-bottom disk (uhpc/air) at `z_ballast`, also when `z_ballast` is at a joint. No face where
  uhpc meets uhpc inside a module, so none at `z_ballast` ≤ `z_lo` (the inner patches close the void).
- Thin shell faces: no joint caps (lateral faces are only split at the module planes);
  `ballast_top` is the full section at `z_ballast` (ballast/shell), split into an annulus
  (ballast/shell) and a disk (ballast/air) when the inner set is cut.
- `brep.bodies`: modular precast `<hull>_UHPC_module_<i>` (solid, one closed shell each; a hollow
  module is a ring); thin shell `<hull>_STEEL_ballast` (solid, full outer section below
  `z_ballast`; absent when `z_ballast` ≤ z_min) and `<hull>_STEEL_shell` (sheet: outer faces
  above `z_ballast`).
- `shells` (signed face lists for F11 and tests): precast `all_outer`, `all_voids{c}` (the fused
  `<hull>_UHPC_all` solid with its cavities); thin shell `layer`, `void` (closed, mass regions only).

**S5 section loop** (F4, F6b): `z` [m, body], `pieces(k)`: `patch`, `u` (iso-u), `curve` (T9 curve
struct, the exact iso-u row in v), `dir` (±1); `pts [n x 3]` CCW on the exact curves (no repeated
point); `area`, `centroid [1x2]`, `I [Ixx Iyy Ixy]` about the origin [m⁴] by Green's theorem (Gauss
on knot spans); `simple`. Body section (F6b): `z`, `module`, `outer`, `inner` (S5 or empty), `solid`.

**S6 body properties** (F6): for each region r and module i: `V(i)` [m³], `S(i,:)` = ∫x dV [m⁴],
`J(:,:,i)` = ∫x xᵀ dV [m⁵], body frame about the origin; `regions.<r>` holds these; `modules(i)`:
`V`, `V_<r>`, `mass`, `rho_eff` = mass/V, `CG_body [1x3]`; `total`: `mass`, `CG_body [1x3]`,
`I_origin [3x3]`, `I_cg [3x3]` (I = ∫ρ(|x|²E − x xᵀ)dV); `quad` (Gauss order used).

**S7 hydrostatics** (F7): `vs`, `draft`, `submersion` (`partial` | `full` | `none`), `V_sub` [m³],
`CB_body`, `CB` (world) [1x3], `S_wet` [m²] (∫|S_u × S_v| over the outer pieces below the
waterline, same Gauss rule), `Aw` [m²], `xF`, `yF` (centre of flotation), `I_wp_xx`, `I_wp_yy`
(about the CF; fixes I12), `I_wp_yy_origin`, `BM_L` = I_wp_yy/V_sub, `KM` (world z) = CB(3) + BM_L,
`waterline` (S5 at z_body = −vs). Waterline at or above z_max: `full`, V_sub = V_hull, CB = hull
centroid, `Aw`, `I_wp_*`, `BM_L` = 0, KM = CB(3), `waterline` empty. At or below z_min: `none`,
V_sub = `S_wet` = `Aw` = 0, CB, `xF`, `yF`, `BM_L`, KM = NaN, `waterline` empty.

**S8 realised design** (`results.stage3`, both modes; `[]` for `preliminary`; replaces
`results.constructability` and `results.steel_data`)
- `mode`, `hull_name`, `status` (`accepted` | `failed`), `reason` ('' when accepted), `escalation`
  (`split` | `fixed_draft` | `spill` | `draft_free`: the last step run), `vs`, `draft` [m] (world).
- `stage2`: `vs`, `rho [N x 1]`, `mass`, `Z_CG` (= `Final3D.CG_total(3)`, world), `GM`
  (`Final3D.GM_L`), `T_heave`, `T_pitch` (`Final3D.periods.heave/.pitch`, coupled).
- `rho` (densities used, by region), `design` (S3), `k_star` (precast: ballast module),
  `V_uhpc_target [N x 1]` (precast split, AGENTS §3 item 4.1). `modules(i)`: `z_lo`, `z_hi` [m,
  body], `t` [m] (NaN: no void), `h_ballast` [m] (from the module bottom, OD6 ii), `V`,
  `V_<region>`, `mass`, `rho_eff`, `rho_stage2`, `rho_floor`, `CG_world [1x3]`.
- `props`: F9 output; `check`: F10 output; `solver`: per escalation step run `exitflag`,
  `iterations`, `fval`, `max_eq_violation`; `fit`: S2r of every inner set; `body`: S4; `step_files`:
  `name`, `path`, `bodies` (F11). `final_props` = `props` + `stage3_status`, `stage3_check`.

## 2. Functions

All take and return the structs of §1; none reads globals or writes files except F11 and F13.

| | Signature | Contract |
|---|---|---|
| F1 | `geo = mwecmass.solid.outer_nurbs(model)`, `outer_nurbs(model, cache, opts)` | S1b from the `.ms2` entity tree. First a deck check on `model.filename` (§8 Scope; `MS2Parser` unchanged): error `UnsupportedEntity` when an entity of the deck that is not in `model.entities` is referenced by the hull (named on the line of an entity the visible surfaces depend on; unreferenced lines are ignored), `UnsupportedMirror` for a mirror plane other than x = 0 or y = 0. Exact path, every entity converted exactly. Points: FramePoint, MirrPoint, AbsBead (on the converted curve, at the parser's parameter map). Curves: BCurve as is, Line degree 1, Arc rational quadratic, BSubCurve by knot insertion, PolyCurve2 joined with C0 knots after exact degree elevation of its pieces to the highest degree, ProjCurve by projecting control points (kept coordinates copied bitwise), EdgeSnake as the boundary row or column of its converted parent surface, taken from the same array (bitwise seam; edge numbers are the parser's, S1: 1 = v0, 2 = u1, 3 = v1, 4 = u0). Surfaces: RevSurf as profile × rational arc (multiples of 90° use exact 0, ±1), RuledSurf as degree 1 in v between its two curves after one common reparameterisation (one knot vector, equal weights), mirrors by flipping control points. `offset_kind` = `rev_z` when the axis end points have equal x and y within T1's 16-ulp bound; the revolution then uses the (x, y) of the axis Line's second end point, every profile control point that is an axis end point by entity (the profile ends at a point or bead that defines the axis) takes that (x, y) bitwise, and control points are axis point + radial vector (a zero radius gives the axis point bitwise); `ruled_parallel` when the cross product of every ruling with the first is bitwise zero. Order of the parameters: when the converted patch is not `z_of_u` but every control column has one z (bitwise) and the weights are separable, z depends on v alone; F1 then exchanges u and v exactly (control points and weights transposed, degrees and knot vectors exchanged; no coordinate changes), sets `swap_uv`, and takes the parser's edge e of that patch as its S1 boundary 5 − e (EdgeSnake, seams); so z depends on one parameter, which F1 orders as u. A patch takes this exact path when its entity converts exactly, it is `z_of_u` (after that exchange) with z(u) monotone, no vertex lies inside one of its rows (§8 Cuts), and every boundary it shares end to end with another exact patch is one converted curve (bitwise, as through EdgeSnake); every other patch takes the general path of §8 (faces fitted through exact points, `exact` false, `fit` filled; exact patches are then split at the band ends). Order of the rules, so that the outcome does not depend on patch names except in the last step: (1) §8 Cuts and Flat regions first; a boundary that holds a neighbour's corner inside it (a T-junction, a boundary shared with only part of another) is resolved only here, in two steps: (1a) a constant-z patch that has a vertex inside one of its rows or a corner inside a row of another patch joins a flat region, and so does every constant-z patch that touches that region at its height (§8 Flat regions); (1b) the vertices are then counted with the end points of seams between patches merged into one region dropped, and a patch that still has a vertex inside one of its rows takes the general path; a seam with a merged flat region is never a conflict. (2) Tie-break, only for a boundary that two exact-convertible patches share end to end (the same two end points) and that is not one curve bitwise: the later of the two in `model.visible_surfs` order takes the general path. (3) The §8 Seams rule for fitted faces. A mirror takes its source's path, so its source and the source's other mirrors go with it; (1)–(3) are repeated until nothing changes. Constant-z intervals: z(u) of a `z_of_u` patch is rational on each knot span, so an interval of u on which z is constant (a horizontal shelf of a profile) is a union of whole spans and ends at existing knot values; F1 splits every exact patch at the ends of each such interval whose height is strictly inside the patch's z range, by knot insertion of that existing value to multiplicity du (none at a `c0_u` row; no root search, no F3b call). Each interval piece is a constant-z S1 piece (`z_range(1)` = `z_range(2)`: F3b never cuts it, F4 skips it, F6 adds 0, as the box's top), or part of a flat region when a vertex lies inside one of its rows (§8). So on every face F1 returns, z(u) = z has one root for every z strictly inside its z range (constant intervals remain only at its lowest or highest z, as the cylinder's end disks; C1: no such interval). A patch boundary that is neither a pole nor shared with another face or a flat region (an open hull, e.g. a hull patch of an entity type the parser skips; the message lists the deck's unparsed lines) raises `HullNotClosed`. Optional inputs `geo = outer_nurbs(model, cache, opts)`: `cache` (T1's) and `opts.t_min` (ε of the fitted faces) are needed only when a patch takes the general path; missing there: error `FitInputMissing`. `opts.max_passes` (default: the kernel's one refinement cap, shared with F2 and set at the T3 checkpoint, §9 item 2): cap reached on a fitted outer face: error `FitNotConverged` listing the failing faces and metrics (M3, `fit.dev_max`). `opts.force_general` (default false) sends every patch through the general path (plan T2b tests it on C1 against the exact path; needs `cache` and `opts.t_min`). C1: all 8 exact, so `outer_nurbs(model)` returns the same `geo` as before. C1 results: §3, C1 oracles. |
| F2 | `[inner, rep] = mwecmass.solid.offset_surface(model, cache, geo, t, z_range, opts)` | S2 at thickness t over `z_range` (body). Nodes from `MS2Parser` + T1 `surface_normals`, offset by d, folds trimmed (`trim_fold.m`), faces split at creases into pieces that meet along u- or v-boundaries, cubic fit (`fit_bspline_surface.m`), knot insertion where M1–M3 fail (judged between the faces as written, §8), knot removal while they hold (AGENTS §5 item 9). Keeps the outer v-construction for `rev_z` (offset profile fitted in u, revolved with the same rational arc) and `ruled_parallel` (both boundary curves offset at the same u nodes, one knot vector, ruled in v), as I3 and §8, when every trim curve of the piece is a parameter line (C1, the cylinder, the box); every other inner piece (`offset_kind` `''`, or trimmed along a curve that is not a parameter line) takes the general path of §8. Inner patches of a `z_of_u` patch are `z_of_u` with monotone z, split at constant-z intervals as in F1. Source points are taken wherever needed so the inner patches cover `z_range`; past an open end (precast: the bottom of the wall module) F5 cuts them with F3b. At a tangent-continuous seam a seam node is offset once, so neighbours share boundary curves bitwise (I2). At a C0 seam (a crease, along u or v) each side is offset along its own normal: on a convex crease the two offsets overlap and are trimmed at their intersection (as a fold); on a concave crease the gap between them is closed by a face of its own, the crease curve offset by d along the fan of normals from one side's to the other's (rule 3: normal distance d everywhere). Where creases meet at a vertex, its offset is the vertex moved by d along every direction of its cone of normals (the directions n with n·w ≤ 0 for every direction w that leaves the vertex out of the solid: one direction at a smooth point, the fan at a concave crease, none at a convex crease or convex vertex); where that cone spans a solid angle (as where every crease at the vertex is concave, e.g. the foot of a re-entrant vertical edge of an L-, T- or cross-shaped column standing on a wider pontoon, where three concave creases meet), the gap between the fan faces of its creases is closed by a face of its own: the part of the sphere of radius d about the vertex over that cone, fitted with z as one parameter as in §8 (it may end in a pole row at its lowest or highest point), its boundaries the end arcs of the adjacent fan faces, shared bitwise (I2); where it overlaps another offset face, the two are trimmed at their intersection as at a fold. `opts.t_min` (eps_fit), `opts.max_passes` (default as F1, §9 item 2); cap reached: error `FitNotConverged` listing the failing patches and metrics. `opts.knots_from` (an S2 set): same pieces and knot vectors, only the control points refitted at the new t (§7 item 4); M1–M3 reported in S2r, not refined. Tests: identity at the same t; properties smooth in t. |
| F2b | `d_close = mwecmass.solid.void_closing_distance(model, cache, geo, z_range)` | Smallest offset distance at which offset layers from opposite sides meet inside `z_range` (not a fold of one layer). C1 neck: 0.10 m. The design bound is t_max = d_close − eps_fit/2: at t = t_max the void closes (d = d_close); an evaluation with d ≥ d_close errors `VoidClosed` (§8). |
| F3 | `[S, Su, Sv] = mwecmass.solid.eval_bspline_surface(surf, u, v)`; `[C, Cs] = …eval_bspline_curve(curve, s)` | Rational or not; u, v column vectors; [n x 3]. |
| F3b | `[lo, hi, u_star] = mwecmass.solid.split_bspline_surface(patch, z)` | `z_of_u` patch, z strictly between z(u0) and z(u1), in either order: u* from z(u*) = z, knot insertion to multiplicity du; `lo`, `hi` keep the parent parameter and share the cut row bitwise. Errors: `ZNotOneParameter` (not `z_of_u`), `ZOutside`, `ZNotMonotonic` (z(u) = z has more than one root). Every face F1 and F2 return is `z_of_u` with one root of z(u) = z for every z strictly inside its z range (exact path by test and by F1's split at constant-z intervals, general path by construction, §8), so `ZNotOneParameter` and `ZNotMonotonic` signal invalid input. |
| F4 | `loop = mwecmass.solid.slice_bspline_surface(patches, z)` | S5 from the iso-u rows; one closed simple loop or error `SectionNotClosed` (as T1). Constant-z patches (z_range(1) = z_range(2)) are skipped. At the height of a flat part inside the hull's range (§8 Scope) the faces ending there from below and from above give two loops; the caller passes the faces of one side (z range below z, or above), as F5, F6b and F7 do. |
| F5 | `body = mwecmass.solid.build_body(geo, design, inner)` | S4. `inner` holds one S2 set per distinct t whose `z_range` covers its modules (§7). Builds every face once; edges shared by index; each flat region (`geo.flat`, `inner.flat`, §8) as one `plane` face (role `outer` or `inner`) whose loops (an outer loop, plus a hole loop where the hull continues through it, as on a shelf) are chains of the lateral-face rows that name it (`[0 j]`), built as caps are. Layout of S3, including `z_ballast` ≤ `z_lo`; faces per mode as in S4. A patch is cut only at planes strictly between z(u0) and z(u1), in either order (constant-z patches never). At every `joint_step` the loop of the larger t must lie inside the loop of the smaller t, tested on the exact loops without tolerance (no intersection of the two loops, and one point of the larger-t loop inside the smaller-t loop); otherwise error `JointNotNested`, which Stage 3 reports as a failed evaluation. A plane at the height of a flat part of the outer surface or of an inner set (§8 Scope; every face that reaches that height ends there, at a band end, a patch end or the cut at the plane) has loops on each side: F4 of the faces below and F4 of the faces above, outer and, where the S4 layout cuts a void there, inner. Each face F5 puts in that plane (cap, `joint_step` or `ballast_top`, as S4 gives them at any other height) covers exactly the area over which the regions on its two sides, read from the loops of each side, are the two regions S4 separates by that face; e.g. a cap between two solid modules covers the intersection of the two outer loops' areas (so on a deck that steps up over part of the section, with a module edge at the step, the area of the loop from above, shared edges included). The rest of each side's area is a flat part, already a face of the outer or inner surface (a constant-z patch or a `flat` region), and it bounds the module or ballast body on the side away from its normal (below it when `normal_z` = +1). Each such face is bounded by arcs of those loops between the points where two loops meet (vertices, §8 Cuts), so it is one plane face per connected part, with an outer loop and hole loops of existing rows, as caps are; an arc two loops share is one row (shared bitwise), taken once; whether an arc lies inside another loop's area is decided by a point-in-loop test of one interior point of the arc on the exact curves (no tolerance: the arc meets no other loop between its ends). No trimmed plane face is needed, and none is written. |
| F6 | `bp = mwecmass.solid.body_properties(body, rho, opts)` | S6 by the divergence theorem with fields (f,0,0): V = ∮x n_x, ∫x = ∮x²/2 n_x, ∫y = ∮xy n_x, ∫z = ∮xz n_x, ∫x² = ∮x³/3 n_x, ∫y² = ∮xy² n_x, ∫z² = ∮xz² n_x, ∫xy = ∮x²y/2 n_x, ∫xz = ∮x²z/2 n_x, ∫yz = ∮xyz n_x. Horizontal faces (caps and constant-z outer patches) contribute exactly 0 (n_x = 0), so only lateral faces are integrated; a region's integral sums faces with `inside` = r minus faces with `outside` = r. Gauss–Legendre of order `opts.n_gauss` on every knot span (exact for polynomial faces, convergent for rational ones). |
| F6b | `sec = mwecmass.solid.body_section(body, z)`, `body_section(body, z, side)` | Body section (S5 block). Each height belongs to the faces and the module above it (half-open [z_lo, z_hi), as SK's stand-in), z_max to those below. `side` (added; default `'above'`, `'below'` the other): at the height of a flat part of the outer surface or of an inner set inside (z_min, z_max) (§8 Scope), F4 of the faces of that side only (z range above z, or below; F4's rule) and `module` the module on that side of a module edge there; at any other height `side` selects in the same way at a module edge or at `z_ballast` (the module and the faces of that side, so the inner loop of that side's t) and has no effect elsewhere. No `SectionNotClosed` at a flat part's height. |
| F7 | `hs = mwecmass.solid.hydrostatics_at_draft(geo, vs, opts)` | S7 on the outer patches (S1, exact or fitted, §8) cut at z_body = −vs (F3b, F6 on the outer pieces, F4 for the waterplane; at a flat part's height, F4 of the faces below); `full` and `none` as in S7, without a cut. |
| F8 | `fl = mwecmass.driver.density_floors(model, cache, geo, edges, t_min, rho_solid, rho_air, solid_modules)` | Per module ρ_min = [ρ_solid·V_solid + ρ_air·V_air]/V with a t_min shell, no ballast (F2, F5, F6); solid modules give ρ_solid. Returns `rho_min`, `V`, `V_solid`, `V_air` [N x 1], `fit` (S2r). |
| F9 | `props = mwecmass.realise.evaluate_realised(bp, hs, design, config)` | `final_props` fields of `export_schema.m` from the exact body. Computed here from the exact body: `vertical_shift` = `design.vs`, `mass_total`, `CG_total` = [0 0 CG_body(3)+vs] (declared symmetric model; x, y kept in `bp`), `Inertia_Tensor` = `I_cg`, `Ixx/Iyy/Izz`, S7 fields, `A_sub` = `hs.S_wet`, `GM_L` = KM − CG_total(3), `periods.heave_uncoupled` = 2π√((M+A33)/K33) and `periods.pitch_uncoupled` = 2π√((I_cg,yy+A55)/K55), Inf when K ≤ 1e-6 (`properties_3d.m` §11; today from `steel_data`), `densities_at_nodes` = `realised_strip_density` = `[bp.modules.rho_eff]'` (layout of x(2:end)), `cross_section` = `config.profile` (as `optim/run.m`), `realised_strip_edges` = `design.edges`, `components(i)` (`density` = `rho_eff`, `z_level` = module mid-height + vs, world, as Final3D), `density_profile_source` = `'realised_partition'`. As `build_realised_properties.m` today: `K_hydro` (K33 = ρ_w g A_w; K55 = M g GM if GM > 0, else 0, `properties_3d.m` §8), A = `mwecmass.bem.interpolate_at_draft(vs, config, CG_total(3))`, `periods.heave`, `.pitch`, `.surge` (Inf, K11 = 0), `coupled_periods`, `participation_factors`, `coupled_modes`, `surge_per_pitch` (`coupled_periods_by_share`), `mass_discrepancy`, `mass_buoyant_force` = ρ_w V_sub, `K_pto` = zeros(3), `K_total` = `K_hydro`, `MassMatrix_CG`/`_Origin` (`build_mass_matrix`), `fill_method` (renamed only if T0b renames it). `realised_strips` leaves the schema (T5): only the realise files T5 and T7 rewrite read it; its data is S8 `modules`. |
| F10 | `check = mwecmass.realise.check_against_stage2(props, stage2, pct, tol_eq, rho_water)` | `metrics(k)` for `Z_CG`, `GM`, `T_heave`, `T_pitch`: `value`, `stage2`, `rel_dev` = abs(X3−X2)/abs(X2), `limit` = pct/100, `pass`; `equalities(k)`: `flotation` (M/(ρ_w V_sub) − 1, as Stage 2) and `GM` (GM/GM2 − 1), `residual`, `tol` = `tol_eq` (the Stage-3 fmincon ConstraintTolerance, 1e-6 today), `pass`; `pass`, `failed` (cellstr), `reason`. |
| F11 | `files = mwecmass.output.step.export_stage3(realised, out_dir)` | Writes with T9 `write_step` into `Output/<type>/step/`: precast `<hull>_UHPC_module_<i>.step`, `<hull>_UHPC_all.step`; thin shell `<hull>_STEEL_ballast.step`, `<hull>_STEEL_shell.step`, `<hull>_STEEL_all.step` (both bodies, junction curve shared). Returns `step_files`. |
| F12 | `data = mwecmass.output.figures.realised_section_data(realised, z_plan, n_z)` | y = 0 elevation (points where F6b loops cross y = 0, found on the exact curves), plan loops at `z_plan`, `z_ballast`, module edges, waterline, `status`. Body frame plus `vs`. F12 calls F6b with two inputs (side `'above'`) at every height, and with `side` `'below'` as well only at the height of a flat part inside (z_min, z_max) (the heights of `geo.flat` and of constant-z pieces there; none on the SK fixtures): a plan loop there is drawn for both sides, and an elevation sample there gives the points of both sides, so a shelf shows as a horizontal step. |
| F13 | `plot_modular_precast(realised, config)`, `plot_steel_solve(realised, config)`, `plot_optimised_cross_section(final_props, config, x_opt, realised)` | Draw F12 data only. |
| F14 | `[results, final_props] = mwecmass.realise.<type>.run(config, x_opt, opt_results)` | Signature unchanged. Sets `results.stage3` (S8) and returns `final_props`; calls F13 and F11 (`config.output.save.stage3.step`) for every status, closest fail included. |

## 3. Ownership and stand-ins

| Item | Implements | Consumers | Stand-in before the join |
|---|---|---|---|
| S1, S1b, F1, F3, F3b, F4 | T2a (general path: T2b, §6) | T3, T4b, T5, T7, T8, T10 | SK fixtures (below) |
| S2, S2r, F2, F2b | T2a (general path: T2b, §6) | T3, T4b, T5, T6, T7 | SK |
| S3–S6, F5, F6, F6b | T3 | T4b, T5, T6, T7, T8, T10 | SK |
| S7, F7 | T3 | T5, T6, T7 | SK |
| F8 | T4b | build_config, `stage2_bounds.m` | — |
| S8, F9, F10, schema of S8 | T5 | T6, T7, T8, T10 | SK `sti_realised`, `sti_stage2` |
| F11 | T10 | F14 | none: call gated off in tests |
| F12, F13 | T8 | F14, `dispatch.m` | none: figure switches off in tests |
| F14 | T5 (precast; T6 adds the escalation inside `solve.m`), T7 (thin shell) | driver | — |

Where the plan's task sections differ, this contract wins:
- Paths: F10 is `+realise/check_against_stage2.m` (both modes, §5 item 11); F9 is the new
  `+realise/evaluate_realised.m`, replacing `build_realised_properties.m`; T5 and T7 write the F14
  calls of F11 and F13, so T8 and T10 do not edit `run.m`.
- Methods: F6 by the divergence theorem on B-rep faces and F7 on the outer patches cut at the
  waterline, not quadrature of sections in z (plan T3); sections (F4, F6b) serve figures and tests.
- Bounds: t_max,i (T6) and `t_max` (T7) are d_close − eps_fit/2 (F2b). Baseline: "updated in the
  same commit" (T4a, T4b) is superseded by §6 (regenerated once, after J2). Switches: one flag,
  `out.save.stage3.step` (F14, §4), not the `out.save.stage3.step_*` switches of the T10 Files line.

**Stand-in kit SK** (one graded step at the start of lane K, merged before any consumer starts):
- `tests/standins/fixtures/` (permanent): `cylinder.ms2` (quarter RevSurf about the z axis of a
  PolyCurve2 of three Lines, bottom disk, side, top disk, mirrored in x and y like C1; sharp rims),
  `box.ms2` (side faces RuledSurf between two vertical Lines: horizontal parallel rulings; bottom
  and top faces RuledSurf between two horizontal Lines of constant z, which F3b never cuts, F4 skips
  and which add 0 to F6: the outer patches close the hull, and F5 keeps one rule, caps only
  strictly inside the z range); `sti_closed_form.m` (V, S, J, sections, hydrostatics for any S3
  design on these, at d: inner cylinder of radius R − d between z0 + d and z1 − d, inner box shifted
  inward by d; `d_close` = min(R, H/2), or half the smallest box dimension); `sti_inner_box.m` (the
  exact S2 set of the box: its six planar ruled faces shifted inward by d, so the void is closed);
  `sti_config.m` (`config` fields `RHO_WATER`, `G`, the one-draft hydro table `hydro_drafts`,
  `hydro_z_cg`, `added_mass_diagonal`, `added_mass_full`, `radiation_damping_full`, `profile`,
  `strip_edges`, `mass_acceptable_pct`; tests put the stand-in S1b into `config.hull_solid`);
  `sti_stage2.m` (Final3D-like, consistent with the closed form); `sti_realised.m` (S8 of a design).
- `tests/standins/+mwecmass/+solid/` and `+realise/`: closed-form versions of F1–F7 including F6b,
  F9, F10 that set or require `geo.analytic` (error `mwecmass:standin:NotAnalytic` otherwise).
- Consumer tests call `addpath(fullfile(root,'tests','standins'),'-end')`: a real function in `src/`
  wins (Octave 8.4 and MATLAB merge package folders, earlier path entries first; verified in
  Octave). The join that merges a producer deletes its stand-in functions; the consumer tests then
  run unchanged on the real code (the join test).
- Fixtures per join test: the cylinder goes through the real F1, F2, F5–F7 (its rims are u-creases
  of convex radius 0 < t, so it also tests fold trimming and crease splitting, plan T2a; rational
  faces, so results are compared with the closed form and printed); the join tests of T4b, T5, T6
  and T7 use it, since their code calls F2. The box serves F5–F7, F11, F12 with `sti_inner_box`
  (its vertical corners are convex C0 v-seams; since §8, general hulls, F2 builds them too, and T2a
  compares its box set with `sti_inner_box`); results polynomial, asserted exact.
- C1 oracles (by construction): `BeadBottom` evaluates to x = 3.9e-16 m, `BeadTop` exactly to
  (0, 1); both ends of `curve6` lie at these beads, so `surface1` (`rev_z`, `z_range` [1.1 −3.25])
  has poles at its keel and neck-top rows. `Edge_For_Dev` = EdgeSnake edge 3 of `surface1` (v1,
  φ = 360°, the profile in y = 1); `curve7`, its projection on y = 0, copies x and z, so `surface2`
  is `ruled_parallel`, and its ridge row u0 (z = 1.1) and keel row u1 have x = 0, shared with its
  X-mirror's. `c0_u` of both: the arc-to-`curve1` joint, z = 1.0; `c0_v` empty; all 8 `exact`;
  v-seams tangent-continuous (`surface1`–`surface2` by construction, mirror seams in symmetry planes).
  For z ∈ [−0.5, 1.0] the section is a stadium of radius 0.10 about (0, ±1), flat sides x = ±0.10
  (curve1's first span has collinear control points), area 0.4 + 0.01π m², inner radius and
  half-width 0.10 − d; z ∈ [1.0, 1.1] is a half cylinder (r 0.10, length 2) with two quarter
  spheres, V = 0.01π + (2/3)π·0.001 m³; neck `d_close` = 0.10 m. Smallest convex principal radius
  0.100 m (top arc, neck hoop), so the offset folds only from d = 0.100 m on (AGENTS rule 3).
- **SK acceptance** (reproduced by the grader): every stand-in S4 (split at `c0_u` rows) passes
  `validate_brep`; the box and cylinder bodies written by `write_step` pass `step_check.py` (closed
  solids, solid count), `occ_volumes` compared with `sti_closed_form` (printed); V, S, J of the
  closed form equal an independent Gauss–Legendre integration in z of the closed-form sections
  (piecewise polynomial, so exact; asserted to rounding); region volumes sum to V_module
  (asserted); every stand-in errors `NotAnalytic` on the C1 `geo`.
- **SK2** (plan; after the general-hulls amendment, before J1): the stand-ins follow every change
  since SK that they mirror: F1 accepts `outer_nurbs(model, cache, opts)` and ignores `cache` and
  `opts`; F6b accepts `body_section(body, z, side)`: `'above'` (the default) is the stand-in's
  half-open rule as it is, `'below'` takes the faces and the module below a module edge or
  `z_ballast` (neither fixture has a flat part inside its z range); S1 entries carry `visible`
  (outer entries: k; S2 pieces: the `visible` of the outer
  entry they are offset from), `fit` (`[]`) and `swap_uv` (false), and S1b and S2 sets carry `flat`
  (empty); the SK field-set assertions list them; F2 returns `sti_inner_box` for the box, so the
  stand-in's `FitNotConverged` branch for the box and the SK test that expects it are retired.

## 4. File ownership map (files more than one task edits; order of edits)

A task edits a shared file only after every earlier task in its row is merged (or, inside the J2
group, accepted) and it has rebased on it; such edits are the task's last commit.

| File | Order |
|---|---|
| `src/+mwecmass/+driver/build_config.m` | T0b → T0d → T4a → T4b → T5 (precast D6 forwards) → T7 (thin-shell D6 forwards) |
| `WEC_User_Input.m` | T0b → T0d → T4a → T5 (precast `t_init`, `max_slope_factor`, `n_z_grid`) → T7 (thin-shell ones) |
| `WEC_Output_Options.m` | T0b → T10 (`out.save.stage3.step`) |
| `src/+mwecmass/+driver/run.m` | nobody (F11 uses `mwecmass.output.output_dir`) |
| `optim/stage2_bounds.m` | T4a → T4b |
| `optim/run.m`, `optim/report_assemble_results.m` | T4a (`run.m`) → T5 (`stage3` replaces `constructability`/`steel_data`) |
| `output/export_schema.m`, `check_export_schema.m`, `export_results.m`, `load_results.m`, `Report.m` | T0b → T5 (S8 for both modes) → T11 (D9 remnants) |
| `docs/RESULT_SCHEMA.md` | T0b → T5 → T12 |
| `docs/RUNTIME_GUIDE.md` | T0b → T4a → T12 |
| `docs/METHODS_ENGINE.md`, `README.md`, `AGENTS.md`, `_graph/CODE_MAP.md` | spec → T0b → T12 |
| `+realise/+modular_precast/*` | T0b → T5 → T6 |
| `+realise/+thin_shell/*`, `+realise/empty_realised_properties.m` | T0b → T7 |
| `+realise/build_realised_properties.m` | T0b → T5 (stops calling it) → T7 (deletes it) |
| `+output/+figures/plot_modular_precast.m`, `plot_steel_solve.m`, `validation/diagnostics/stage_animations.m` | T0b → T8 |
| `output/dispatch.m` | T4a (only if the `Monotonic density` line goes with `c_mono`) → T8 (F13 call) → T10 (`validate_save_flags`: `out.save.stage3.step`) |
| `tools/baseline_run.m` | T0 → T5 (summary reads `results.stage3`) |
| `validation/diagnostics/uhpc_mass_balance.m` | T0b → T11 |
| `tests/run_tests.m`, `.gitignore` | T0 → T0d |
| `tests/baseline/*.json` | T0 → J2 run step |
| `tests/standins/*` | SK → SK2 → J1 (deletes the stand-ins of F1–F7, F6b) → J2 (deletes the F9, F10 stand-ins); consumers do not edit it |
| `tests/standin_kit/*` | SK → SK2 (retires the box `FitNotConverged` assertion, extends the field-set assertions) → J1 (removes or adapts the SK tests of the stand-ins it deletes) |
| `MS2Parser.m` | T0c only. `+solid/*`: T2a, T2b and T3 own disjoint files (T2b: `fit_z_faces.m` and its `fit_z_*.m` helpers, the general-path fixtures in `tests/solid/fixtures/` and its tests; T3 adds `tests/solid/fixtures/stepped_box.ms2`; J1 adds `tests/solid/test_join_general.m`). `+output/+step/*`: T9 → T10 (new files; `write_step.m` changes only if a C1 body fails `step_check`, by T3). |

## 5. Invariants every implementation tests

- **I1 Volume closure (rule 11)**, to a stated multiple of machine epsilon justified in the test:
  per module V_uhpc + V_air = V_module (precast), V_ballast + V_shell + V_air = V_module (thin
  shell); Σ_i V_module,i = V_hull, V_hull integrated on the unsplit outer patches with the cut
  parameters u* as extra breakpoints (same (u, v) samples). With `z_ballast` inside a module, at a
  module edge, spilled, and below the inner `z_lo`.
- **I2 Shared edges:** every closed shell passes `validate_brep`; the boundary rows of the two faces
  at an edge equal its curve bitwise (`isequal`: control points, knots, weights), on u- and v-seams,
  including C1's keel seam of `surface2` and its X-mirror (§3, C1 oracles); on a fitted hull, no
  vertex lies inside a face row (§8 Cuts), so the seam fields pair up (face k boundary b names
  [j c] exactly when face j boundary c names [k b]).
- **I3 Exact z-splits** where z depends on u only. RevSurf about a vertical axis: rotation keeps z,
  so each control row has one z and w_ij = a_i b_j, hence z(u,v) = Σ N_i a_i z_i / Σ N_i a_i. C1
  RuledSurf: `curve7` keeps z at both ends of every ruling. A plane z = z_k meets the patch in the
  iso-u row u*; knot insertion cuts it without moving the surface (test: parent and pieces agree at
  the same (u, v) to rounding; cap vertices within T9's 1e-7 m). Offsets keep it for `rev_z` and
  `ruled_parallel` only (normal in the meridian plane; normal constant along a ruling). Every other
  face, outer or inner, is fitted with z as one parameter (§8), so its cut is the same knot insertion.
- **I4 Units:** METRE everywhere; `step_check` reports `declared_length_unit = METRE`. **I5 Names:**
  new code uses only the names of §0 (`grep` for the old names in new files is empty).
- **I6 Frames:** CG_total(3) = CG_body(3) + vs; KM, CB world. **I7 Orientation:** outer normals
  out of the hull, inner normals into the void; every region volume positive.
- **I8 C1 symmetry:** CG x = y = 0 and products of inertia 0 to machine precision.
- **I9 Exactness of F1:** F4 sections of the NURBS and T1 `outer_rows` at the same z: the
  difference is printed per height and asserted within T1's documented in-plane accuracy at that
  height (about sqrt(16 ulp/|z_ss|), 5.5e-8 m at the C1 shoulder z = −1), not to rounding. This holds
  for exact patches; a fitted face instead has `fit.dev_max` ≤ ε/4 (§8), asserted and printed per face.

## 6. Lanes and joins

```
J0   merge T0, T1, T9, this contract
K    SK ─► T2a ─► T3 ─────────────────────► J1 (T2a+T2b+T3)  ── owner checkpoint after T3
            T2a ─► T2b (∥ T3; both from T2a's branch, disjoint files) ─┘
     SK2 (stand-ins and SK's tests only): right after the general-hulls amendment merges, before J1
N    T0b ─► merged at once (all later lanes start from it)        T0c (MS2Parser.m) in parallel
P    [T0b, SK, T0c merged] ─► T0d ─► T4a ─► T4b (after SK2; stand-ins; merges after J1)
U    [T0b, SK merged] ─► T5 ─► T6
S    [T0b, SK merged] ─► T7
O    [T0b, SK merged] ─► T8 ∥ T10
J2   after J1 and lane P: the group T5, T6, T7, T8, T10 merges as one chain in that order
     ─► first whole-pipeline runs ─► owner checkpoints after T6 and after T10
G    after J2: T11 ─► T12 ─► T13
```
- T2a (exact path) is the critical path to T3; T2b (general path) runs beside T3. J1 merges T2a,
  T2b and T3, in that order, deletes the stand-ins of F1–F7, F6b, and adds the join test of the
  general path through the real F5–F7 (plan J1).
- SK2 changes are additive for lanes U, S and O (new fields at their defaults; a new optional F1
  form they do not call; the optional F6b `side`, which T8's F12 passes only at a flat part's height
  inside the hull's range, which no stand-in fixture has; F2 on the box, which they do not call
  either, since the box serves them
  through `sti_inner_box`; F2 faces at concave creases and vertices, and the F4, F5 rules for flat
  parts, which no stand-in fixture has), so those
  lanes need not wait for it. T4b calls F1 with three inputs (S1b for any deck, §1) and starts after SK2.
- Tasks are graded on component tests (stand-ins, C1 pieces, real kernel after J1). Whole-pipeline
  items are checked after J2 (owner): the T0 baseline on the commit after T0d (identity of T0b, T0c,
  T0d); C1 runs of both modes on the J2 head for T4a, T4b, T5, T6, T7, T10; then the baseline is
  regenerated once (`tools/write_octave_baseline.m`), changes printed per task. Failures are fixed
  on a fix branch; inside J2 the integration branch is pushed after the last merge and suite.

## 7. Stage-3 evaluation per solver iteration, and what must be cached

Per evaluation of (t_k*, t_i, `z_ballast`, vs) or (t, `z_ballast`, vs): S2 sets per distinct t (F2,
the expensive step), F5 (knot insertion at `z_ballast`), F6 on the changed faces, F7 (vs only), F9
(BEM interpolation, 3 × 3 eigenproblem), F10. Needs:
1. S1b, the outer pieces cut at the module edges and their F6 integrals: once per deck and edges.
2. S2 sets keyed by the bitwise value of t only, each fitted once over the full hollow range
   (precast: z_min to the bottom of the solid wall module; thin shell: the whole hull): one set
   serves every module with that t (equal t shares joint edges; `z_ballast` moves without
   refitting), and a finite-difference step on one t never refits another.
3. F7 keyed by vs: one call while the draft is fixed (`split`, `fixed_draft`, `spill`).
4. Smoothness: properties differentiable in t and `z_ballast` within one solve, so each S2 set keeps
   its knot vectors during a solve (F2 `opts.knots_from`); the adaptive fit runs on the final t, and
   a changed knot structure restarts the solve there. `z_ballast` enters only by knot insertion.

## 8. Decisions taken here (simplest option that meets the rules)

- One B-rep (S4) for mass, sections, figures and STEP: one geometry (§5 item 1); caps add 0 to F6.
- Inner fit keeps the outer v-construction only where it is the exact normal offset (I3): exact
  z-splits, no trimming curves in the T9 writer. `rev_z` within T1's 16-ulp bound and the F1 axis
  snapping: no new tolerance, and C1's poles and keel seam come out bitwise.
- One S2 set per distinct t (§7 item 2). `A_sub` = F7 `S_wet` on the outer pieces (one kernel).
- No margin below t_max (agent's decision, owner informed): solvers use t_max = d_close − eps_fit/2;
  an evaluation with d ≥ d_close returns error `VoidClosed` and counts as a failed evaluation.
- `results.stage3`, F9, F10 shared by both modes, schema owned by T5 (one struct, one owner;
  OD11); the thin-shell process stays its own (rule 7). `<hull>_STEEL_all.step` mirrors `_UHPC_all`.
- **General hulls** (owner's decision of 2026-10-07 on former §9 item 1: "No it needs to be
  generalized."; AGENTS §3 item 36). It follows AGENTS §5 item 9: exact NURBS where the entity
  allows, otherwise fitted with the same metrics.
  - *Scope (one line).* Any entity type `MS2Parser` evaluates and any parametrisation, for hulls
    whose horizontal section is one closed loop at every height (one body, no holes); a section of
    several loops raises `mwecmass:solid:SectionNotClosed` (T1, F4).
    Details: every height in (z_min, z_max), as T1 `outer_rows` already assumes, except the height
    of a flat part inside that range (a shelf: constant-z pieces or a flat region), where the hull
    is horizontal over an area that the loops just below and just above bound (F4, F5). Parser boundary: `MS2Parser` (unchanged, T0c)
    evaluates FramePoint, MirrPoint, AbsBead, AbsRing, BCurve, Conic, CopyCurve, Line, BSubCurve,
    Arc, PolyCurve2, ProjCurve, EdgeSnake, BSubSnake, RuledSurf, RevSurf, BLoftSurf, DevSurf and
    MirrSurf, and mirrors only about x = 0 or y = 0; it skips every other entity type without a
    message (including `RealList`, `Variable` and `Pathname`, on purpose) and reads any other
    mirror plane as y = 0 with only a warning. F1 therefore checks the deck first (re-reading
    `model.filename` with the parser's line rules) and stops only a hull the parser would misread:
    `mwecmass:solid:UnsupportedEntity` when a deck entity that is not in `model.entities` (no parser
    case, or its line failed to parse) is referenced by the hull, i.e. named as a token on the line
    of an entity in the closure of the visible surfaces (each visible surface, a mirror's source,
    and recursively every entity named on the line of an entity already in the closure); a line
    nothing in the closure names (display or analysis entities, skipped types used by nothing) is
    ignored. A hull patch whose own entity the parser skips is in no closure; the hull it leaves
    open raises `mwecmass:solid:HullNotClosed` (F1: a patch boundary with no neighbour that is
    not a pole). `mwecmass:solid:UnsupportedMirror` (T1's identifier) is raised for a MirrSurf
    plane or a `Symmetry:` entry other than x = 0 or y = 0, such as a hull mirrored top to bottom.
  - *Two paths per patch.* Exact path (F1; conditions there): the entity converts
    exactly and z depends on one parameter (F1 orders it as u, exchanging u and v exactly where z
    depends on v alone), monotone (all of C1). General path: every other patch
    is replaced by faces fitted through exact points of the parametric definition, the T1
    `outer_rows` sections (points on the exact surface with their patch and (u, v)) at heights
    placed adaptively in z (AGENTS §5 item 9.1), with T1 `surface_normals`. A face covers one
    height band, in which its part of every section is one arc between two points of its
    v-boundaries (seams, creases or cuts). The band ends are one set of heights for the whole
    hull: every patch corner, every z-extreme of a patch boundary or crease curve, every row where
    the surface is horizontal (where z may turn back), every single highest or lowest point inside
    a patch (the faces end there in a pole row, S1 `pole`; the point solves z_u = z_v = 0 on the
    exact surface), and z_min, z_max. Every face, fitted or exact, is split at every band end
    strictly inside its z range (an exact patch by F3b's knot insertion, which keeps it exact; the
    pieces replace it in S1), so every seam, crease or cut between two faces spans one band. A band
    end at the height of a constant-z interval of an exact patch (where the surface is horizontal
    over an interval, e.g. the step of a stepped profile revolved about a vertical axis) is never
    strictly inside a piece: F1 has already split the patch at the existing knots that bound the
    interval (F1, constant-z intervals; no F3b call), and the interval piece is a constant-z S1
    piece, or part of a flat region if a vertex lies inside one of its rows. In a
    face u runs with z and v along the arc: every control row has one z (bitwise) and the row
    heights are monotone in u, strictly except for a first or last pair of rows that meets a
    horizontal tangent, so the face is `z_of_u` with monotone z by construction. Every face is an
    untrimmed patch, and every horizontal cut (module edges, `z_ballast`, waterline) is a parameter
    line found by F3b's knot insertion, so the T9 writer, F5–F7 and the SK stand-ins need no
    trimming curves. Fitted faces are non-rational cubics except where a seam imposes a
    neighbour's structure (Seams, below). Inner pieces whose structure is not kept (F2) are fitted
    the same way through their offset nodes, and so are the fan faces of concave creases that keep
    no such structure (the fan of a u-crease of a `rev_z` piece is a revolved arc and keeps it) and the
    sphere faces of vertices whose cone of normals spans a solid angle (F2), which close the inner
    surface where concave creases meet; every rule of this item applies to an S2 set on its own.
  - *Cuts: every face boundary has one neighbour (no T-junctions).* A vertex is a point of a
    band-end row where a v-boundary of a face (a seam, crease or cut) meets or ends, where a face's
    row collapses (a pole, or the narrow end of a face between two curves that meet), where a
    row turns back on itself (an end of a ridge or trough line inside a patch, such as the top line
    of a capsule drawn as one patch), or, at the height of a flat part, where the loop from below and
    the loop from above meet (cross, touch, or begin or end an arc they share; F5 bounds the faces
    in that plane by the arcs between these points). A vertex is a point in space: it is a vertex of every face row
    through it, on both sides of the band end and on both sides of a ridge. Examples in the T2b
    fixture `tilted_revolution`: the highest points of the N, S and Q rims (φ = 180°, where no seam
    runs; each rim's two crease branches end there), and the disk centres T and K, where the self-seams
    of the top and bottom disks end. Faces are cut until no vertex lies inside a row, in two steps.
    1. *Closed-loop regions.* Where a patch's part of every section of a band is a closed loop with
       no seam or crease point (a smooth dome, a crowned deck, a revolution about a non-vertical
       axis whose self-seam is off the loop), the region is cut first along one z-monotone curve on
       the exact surface: the steepest-ascent (or descent) line of z (direction g⁻¹∇z in the patch's
       (u, v), g the first fundamental form). Consecutive closed-loop regions (across band ends with
       no vertex) form a chain cut by one curve, across a row where the surface is horizontal from
       the point where the line arrives. Its start: the vertex of smallest (u, v) on the chain's
       lowest row, else its highest row, at which the surface is not horizontal; if there is none,
       the point of smallest u (then v) of the chain's mid-height section. The cut is built once and
       used bitwise as both the v0 and the v1 boundary of the face (a self-seam, as on a 360°
       revolution: `seam_v0` = [k 3], `seam_v1` = [k 1]; F5 lists its edge once in each direction in
       the face's loop, which T9 `validate_brep` accepts). Where the cut ends on a row at a point
       that is not a vertex, that point becomes one.
    2. *Every other vertex.* A face with vertices inside a row is cut from each of them to a vertex
       of its opposite row: with b_1 … b_m and t_1 … t_n the vertices inside its lower and upper rows
       in the order of v, b_i is joined to t_i for i ≤ min(m, n), and each further vertex to the
       end of the opposite row on the face's v1 boundary. The node of a cut at height z is the
       point of the face's section arc (exact, from T1) at the fraction of the arc's length from
       its v0 end interpolated linearly in z between the fractions of the cut's two ends, so cuts
       are z-monotone, run through exact surface points and do not cross; a cut that ends at a row
       end gives a face with that row collapsed (S1 `pole`), and no face gains a second one.
    Step 2 ends every cut at an existing vertex or a collapsed row and creates no vertex; step 1
    creates at most one at each end of a chain, which step 2 then takes in. Afterwards the rows on
    the two sides of every band end, and of every ridge, run between the same consecutive vertices,
    so every face boundary has exactly one neighbour: the S1 seam fields are unchanged (one
    `[patch boundary]` per boundary) and every shared boundary is built once and taken by both faces
    bitwise (Seams, I2). An exact patch or a structured inner piece with a vertex inside one of its
    rows would need a cut that is not a parameter line: it takes the general path instead, a
    constant-z one by joining a flat region (step 1 of F1's order of the rules, before the
    tie-break; repeated until nothing changes). This is the only rule for a boundary shared with
    only part of another (a T-junction): which patch leaves the exact path follows from where the
    vertex lies, never from the patch names.
  - *Flat regions.* A constant-z region on the general path (a flat deck, a flat bottom, a horizontal
    shelf) is not fitted. One flat region is one connected constant-z area at one height, merged
    across patch seams and mirror planes: the source and mirrors of a flat deck (a mirror takes its
    source's path), a deck drawn as two patches, and every exact constant-z patch or piece that
    touches the area at that height all belong to the one region, which lists them in `visible`.
    It is an entry of S1b `flat` (S2 `flat` for the inner surface) and has no boundary of its own
    with another flat region or constant-z patch: only rows of lateral faces bound it, and each such
    row names it (`seam_*` = `[0 j]`); the seams between the patches it merges are no face
    boundaries, so their end points are vertices (Cuts) only where another rule makes them so.
    Regions are merged before the vertices are counted (F1 order of the rules, steps 1a, 1b): a
    corner of a merged patch on a lateral face's row (e.g. the corners of `cross_column`'s pontoon-top
    rectangles on the walls' top rows, plan T2b) puts no vertex there, and a seam between a lateral
    face and a merged region is never a seam conflict. F5 builds it as one T9 `plane` face at its
    height whose loops
    (an outer loop, plus a hole loop where the hull continues through it, as on a shelf) are chains
    of those rows, as it builds caps, so a vertex anywhere on its boundary is an edge end; F4 skips
    it and F6 adds 0 (n_x = 0). An exact constant-z patch with a vertex inside a row, or with a
    corner inside a row of another patch, takes the general path and joins such a region (F1 step
    1a); the all-exact box keeps its top and bottom as S1 patches.
  - *Mirrors.* A mirror takes its source's path. The mirror of a fitted face is never fitted on
    its own: it is the source's fitted face with the control points flipped (`flips`), knots and
    weights copied, `outward` reversed, as for an exact mirror. A boundary of a fitted source face
    that it shares with one of its own mirrors (T1 `outer_rows` reports those seam points with
    `patch2` = that mirror) lies in the mirror plane: it is fitted with its control points set to
    exactly 0 in the flipped coordinate, so the boundary and its flip are bitwise equal (I2;
    −0 = +0) and the written hull is exactly symmetric; `fit.dev_max` includes this constraint.
    Since a mirror face is never cut on its own, the vertices of a mirror face's rows (those its own
    neighbours put there, which need not be the mirror images of the source's), flipped back, are
    vertices of its source's rows: the source is cut (Cuts) at the union of both sets, so no vertex
    lies inside a row of the source or of any of its mirrors.
  - *Seams (I2).* A boundary shared by two faces is built once and taken by both bitwise. A fitted
    face that meets an exact patch or a structured inner piece takes that boundary curve as it is
    (degree, knots, weights, control points) in that direction and fits only its other control
    points; knots its fit needs are inserted in the neighbour as well, which stays exact. A seam
    between two fitted faces is fitted once. A fitted face that would have to take two different
    structures in one direction from two neighbours (different control-row heights or weights after
    mutual knot insertion) cannot; the later of those neighbours in `model.visible_surfs` order then
    takes the general path as well, repeated until no conflict is left (C1: no fitted face).
  - *Metrics on the written faces.* M1–M3 are judged between the faces actually written, outer as
    written (exact or fitted) and inner as fitted, and F6 and F7 integrate those faces (AGENTS §5
    item 9.7). A fitted outer face is judged by M3 (no self-intersection; every horizontal slice of
    the written outer surface is one simple closed loop, one per side at a flat part's height) and by `fit.dev_max` ≤ ε/4; F1 refines it
    as F2 refines inner faces (`opts.max_passes`; cap reached: `FitNotConverged`). Derivation of
    ε/4, from the owner's band ε = 0.01·t_min alone (rule 5): at a check point let e_o be the normal
    deviation of the written outer face from the exact surface and e_i that of the inner face from
    the exact offset at d = t + ε/2, both positive toward the other face. To first order
    t_local = d − e_o − e_i, so M2 (t ≤ t_local ≤ t + ε) holds for every sign of the errors if and
    only if |e_o| + |e_i| ≤ ε/2, and M1 follows from M2 because t ≥ t_min. The outer faces are
    fitted once per hull, before any t is known, and must leave every inner fit the same budget
    whatever the signs; the two errors enter t_local alike, so each gets half: |e_o| ≤ ε/4. The
    inner fit is judged by M1–M2 on the written faces, which leaves it ε/2 − `dev_max` of its outer
    face: at least ε/4, and ε/2 on an exact patch, as before. No other number enters.
  - *Errors.* `NotExact` (F1) and `OffsetNotStructured` (F2) are removed: those patches and pieces
    take the general path. `ZNotOneParameter` and `ZNotMonotonic` stay in F3b as invalid-input
    guards. New: `FitInputMissing` and `HullNotClosed` (F1). Unsupported, with a named error: a
    section of several loops (`SectionNotClosed`), an entity type `MS2Parser` does not evaluate
    referenced by the hull (`UnsupportedEntity`, F1), and a mirror plane other than x = 0 or y = 0
    (`UnsupportedMirror`, F1 and T1). A plane at a flat part's height needs no error of its own: F5
    builds its faces from the arcs of the loops of both sides, shared edges included. Invalid input:
    an open hull (`HullNotClosed`, F1).

## 9. Open for the owner, and decisions due

1. Decided by the owner on 2026-10-07 (general hulls): moved to §8.
2. Due at the T3 checkpoint: `opts.max_passes` and `opts.n_gauss` (limits, not gates), set from the
   pass counts and convergence tables T2a, T2b and T3 measure and present there.
