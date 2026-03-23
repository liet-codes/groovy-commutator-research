# PLAN.md — Current Status and Next Steps

*Last updated: 2026-03-22 (Cycle 2)*

## Completed
- [x] Core library: ECA engine, operators, measures, affine decomposition
- [x] Experiment 1: Baseline characterization (all 256 rules)
- [x] Experiment 2: Commutator zoo (8 variants × 20 rules)
- [x] Experiment 3: Class IV hunt (extended runs)
- [x] Experiment 4: Affine × commutator analysis
- [x] Paper draft v1: definitions, results, conjectures
- [x] Experiment 5: 88 equivalence class lookup table — confirms all conjectures across 88 classes
- [x] Paper Appendix A: equivalence table with analysis
- [x] Experiment 6: Spectral analysis of G time series — 1/f signatures found for Class IV
- [x] Experiment 7: Glider detection via G peaks — 276 tracks in Rule 110
- [x] Experiment 8: Finite-size scaling (N=51–501) — Rule 54 not converged, Rule 110 stable
- [x] Paper updated with all exp6-8 results (sections 4.5-4.7, updated future directions)

## Key Findings (Cumulative)
1. Class IV rules have highest mean G entropy (0.892) with tightest std (0.079)
2. Class IV has longest correlation lengths (1.63 mean); Rule 54 peaks at 3.50
3. No Class IV rule is affine (0% vs 27% of Class III)
4. Nested commutator [E,[E,D]] drops significantly for Rule 110 (0.973 → 0.713)
5. Multi-step [Eⁿ,D] grows for Class IV but saturates for Class III
6. Conjecture 1 confirmed: all affine rules have G=0 across all 88 equivalence classes
7. **NEW (Cycle 2)**: Class IV rules show positive spectral β exponents (0.25–0.67), suggesting 1/f-like noise; Class III β≈0 (white noise)
8. **NEW (Cycle 2)**: Rule 54 has strongest 1/f signal (β=0.67) — most structured Class IV rule
9. **NEW (Cycle 2)**: 276 glider tracks detected in Rule 110 via G peaks, speeds 0.06–2.10 cells/step
10. **NEW (Cycle 2)**: Rule 54 entropy grows with system size (0.73→0.84 at N=501), not yet converged; Rule 110 converges by N=51

## Immediate Next Steps
1. **Longer spectral runs**: T=10000 or T=50000 to get cleaner β estimates with higher R²
2. **Welch/DFA spectral methods**: Replace raw FFT with Welch's method or detrended fluctuation analysis
3. **Improved glider detection**: Speed-dependent thresholds, cross-correlation tracking
4. **Larger scaling runs**: N=1001, N=2001 for Rule 54 to find its thermodynamic limit
5. **Validate glider tracks**: Compare detected Rule 110 tracks with known glider catalogs

## Medium-Term
6. **Machine learning classifier**: Use commutator features (entropy, β, corr length) for Wolfram class prediction
7. **2D CA**: Extend commutator analysis to 2-state 2D totalistic rules
8. **Algebraic proof**: Prove Conjecture 1 formally (commutator vanishing for affine rules over GF(2))
9. **Connection to universality**: Investigate if nested commutator drop predicts computational universality

## Open Questions
- Why does Rule 54 show growing entropy with system size while Rule 110 converges? Does Rule 54 eventually approach Class III entropy?
- Why does the nested commutator drop for Rule 110 but not Rule 54? Is there sub-structure within Class IV?
- Can the spectral β exponent serve as a standalone Class IV discriminator?
- The 728 "tracks" in Rule 30 — can we develop a false-positive rejection criterion to distinguish genuine gliders from noise?
- Connection to Brooklyn's SAT phase transition: is the commutator measuring constraint density?
- Does the 1/9 universal tiling ratio (Patrick's result) show up in commutator statistics?

## Heartbeat Tasks
When running on heartbeat cycles:
1. Check if any experiments need re-running with updated parameters
2. Try improved spectral methods (Welch, DFA)
3. Run larger scaling experiments (N=1001+) for Rule 54
4. Try new commutator variants suggested in RESEARCH.md
5. Update paper with any new results
6. Always update this PLAN.md when done
