# C1_graph graph report

The graph maps `MASS_SUITE_C1_COMPLETE_R2026a.zip` exactly (197 archive files).

## Counts

- Nodes: 540 (336 functions, 4 classes, 3 runtime modes).
- Edges: 936 (calls=506, contains=340, declares=4, produces=31, references=42, references_handle=7, uses=6).
- MATLAB files: 158; Markdown files: 6; non-code artifacts: 33.

## Scope and limits

The inventory is generated from a clean extraction of the target ZIP. No checkout, Git metadata, legacy tree, tests, caches, generator, or hidden file is consulted. The MATLAB language has no native AST in Graphify, so function declarations and calls are represented by a conservative package-local regular-expression scan; dynamic dispatch, handles assembled at runtime, and ambiguous bare names may be absent. `references_handle` edges record `@name` function-handle evidence only; they are deliberately not `calls` edges.

All paths are archive-relative and all edge endpoints resolve to nodes. HTML is offline-capable and contains no external scripts, stylesheets, or URLs.

## C1 payload relationships

Exactly two data inputs and 31 generated outputs are represented as file nodes. They are grouped under the three runtime modes `preliminary`, `thin_shell`, and `modular_precast`; each mode is linked to both inputs and each mode produces its assigned result/figure/log outputs.
