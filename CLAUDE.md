# CLAUDE.md — Groovy Commutator Research

## Project Overview
Research codebase investigating the Groovy Commutator as a characterization tool for Wolfram Class IV cellular automata. The goal is to find a computational definition of Class IV using commutator structure.

## Project Structure
```
src/           - Core library modules
  eca.py       - Elementary Cellular Automata engine (all 256 rules)
  operators.py - Commutator operators (D, E, I, G, and variants)
  measures.py  - Quantitative measures (entropy, correlation, compressibility)
  affine.py    - Affine decomposition (Brooklyn Rose approach)
experiments/   - Runnable experiment scripts
  exp1_baseline.py        - All 256 rules baseline characterization
  exp2_commutator_zoo.py  - Compare commutator variants
  exp3_class4_hunt.py     - Focused Class IV analysis
  exp4_affine.py          - Affine decomposition × commutator
  exp5_equivalence_table.py - 88 equivalence class lookup table
results/       - Generated CSVs, plots, and tables
paper/         - Formal paper draft (markdown)
PLAN.md        - Living document: current status, next steps, open questions
```

## Coding Practices
- Python 3.10+, numpy for array ops, matplotlib for plots, scipy for stats
- All experiments output to `results/` as CSV + plots
- Use the virtual environment: `source .venv/bin/activate`
- Run experiments from repo root: `python -m experiments.exp1_baseline`
- Keep functions pure where possible; side effects (file I/O) at top level only
- Every experiment script should be independently runnable
- Commit results alongside code (this is a research repo, reproducibility matters)

## Critical Directive
**Always keep PLAN.md up to date.** Before starting work, read PLAN.md. After completing work, update PLAN.md with:
- What was done
- What results were found
- What the next steps are
- Any new questions or hypotheses

This allows any agent (human or AI) to pick up work from where it left off.

## Key Concepts
- **Groovy Commutator**: G(S) = D(E(S)) ⊕ E(D(S)) — measures non-commutativity of spatial differentiation and temporal evolution
- **Class IV hypothesis**: Class IV = nonzero non-affine perturbation AND structured (non-maximal) commutator entropy AND elevated spatial correlation
- **88 equivalence classes**: ECA rules under left-right reflection and black-white complementation
- **Affine decomposition**: Each rule = closest XOR-linear function + nonlinear perturbation residual

## References
- RESEARCH.md — Full research plan and background
- paper/groovy_commutator.md — Current paper draft
- results/equivalence_table.md — Citable lookup table for all 88 equivalence classes
