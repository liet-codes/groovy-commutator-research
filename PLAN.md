# PLAN.md — Current Status and Next Steps

*Last updated: 2026-03-22 (Cycle 3 — Post D operator correction)*

## Critical Change: D Operator Corrected
The differentiation operator D was fixed from spatial adjacency XOR to D(S) = S ⊕ φ(S) (which cells are about to flip). **All experiments re-run with corrected operator.** All results below use the corrected definition.

## Completed
- [x] Core library: ECA engine, operators, measures, affine decomposition
- [x] **D operator fix**: D(S) = S ⊕ φ(S) in `src/operators.py` (commit 695c427)
- [x] Experiment 1: Baseline characterization (all 256 rules) — re-run
- [x] Experiment 2: Commutator zoo (7 variants × 20 rules) — re-run, [E,I] removed
- [x] Experiment 3: Class IV hunt (extended runs) — re-run
- [x] Experiment 4: Affine × commutator analysis — re-run
- [x] Experiment 5: 88 equivalence class lookup table — re-run
- [x] Experiment 6: Spectral analysis of G time series — re-run
- [x] Experiment 7: Glider detection via G peaks — re-run
- [x] Experiment 8: Finite-size scaling (N=51–501) — re-run
- [x] Paper updated with all corrected results, operator correction note, revised conjectures

## Key Findings (Corrected D Operator)
1. Class IV rules have mean G entropy 0.842±0.217; Class II dropped to 0.368±0.361 (many now G=0)
2. 22 of 65 Class II equivalence class representatives now have G=0 (up from ~0 with old operator)
3. Correlation lengths collapsed to ~1.0 across all classes — no longer discriminates Class IV
4. No Class IV rule is affine (0% vs 27% of Class III); Conjecture 1 confirmed (affine ⟹ G=0)
5. Nested commutator [E,[E,D]] drop for Rule 110: 0.938→0.910 (less dramatic than old 0.973→0.713)
6. Some Class II rules have surprisingly high G entropy: Rule 24 (0.974), Rule 35 (0.992)
7. Rule 110 spectral β=0.696 (strongest 1/f signal); Rule 54 β=-0.034 (nearly white noise!)
8. Rule 106 β=0.396, Class IV mean β≈0.352
9. Glider tracks: Rule 110: 564, Rule 54: 394, Rule 30: 770
10. Rule 54 entropy grows with system size (0.735→0.837); Rule 110 stable at ~0.939

## Revised Conjectures
- **Conjecture 1** (confirmed): All affine rules have G≡0
- **Conjecture 2**: Class IV ⟹ G entropy in (0.3, 1.0) with nonzero perturbation weight (threshold lowered from 0.5)
- **Conjecture 3**: Nested commutator [E,[E,D]] entropy < [E,D] entropy for Class IV (still holds)
- **Conjecture 4**: Class IV spectral β > 0.2 (holds for Rules 110, 106; fails for Rule 54)
- **Conjecture 5**: Glider-bearing rules have ≥100 detected G-peak tracks (holds for all tested rules, but high false positives in Class III Rule 30)

## Key Changes from Old Results
| Metric | Old (spatial XOR) | New (D = S⊕φ(S)) |
|---|---|---|
| Class IV G entropy | 0.892±0.079 | 0.842±0.217 |
| Class II G entropy | 0.647 | 0.368±0.361 |
| Rule 4 G entropy | 0.803 | 0.000 |
| Correlation lengths | Class IV elevated (1.63 mean) | ~1.0 across all classes |
| Rule 54 spectral β | 0.67 | -0.034 |
| Rule 110 spectral β | ~0.25 | 0.696 |
| Rule 110 glider tracks | 276 | 564 |

## Immediate Next Steps
1. **Resolve Class II contamination**: Investigate why Rules 24, 35 have G entropy >0.97 (higher than some Class IV rules). May need additional features beyond entropy alone.
2. **Longer spectral runs**: T=10000+ to get cleaner β estimates; current β for Rule 54 is suspiciously low.
3. **Welch/DFA spectral methods**: Replace raw FFT with more robust methods.
4. **Correlation length alternatives**: Since spatial corr length collapsed to ~1.0, try block entropy correlation or mutual information at various lags.
5. **Larger scaling runs**: N=1001+ for Rule 54 to find thermodynamic limit.

## Medium-Term
6. **Multi-feature classifier**: Use commutator entropy + spectral β + perturbation weight + nested commutator ratio for Wolfram class prediction
7. **2D CA**: Extend commutator analysis to 2-state 2D totalistic rules
8. **Algebraic proof**: Prove Conjecture 1 formally (commutator vanishing for affine rules over GF(2))
9. **Rule 54 anomaly**: Understand why Rule 54 shows white-noise spectrum despite being Class IV

## Open Questions
- Why do some Class II rules (24, 35) have very high G entropy? Are they misclassified, or does G entropy alone not suffice?
- Why did correlation lengths collapse to ~1.0 with the corrected D? The old operator created artificial spatial structure.
- Rule 54's spectral β ≈ 0 is surprising — is this a finite-size/finite-time artifact, or genuinely white noise in its commutator?
- Can we combine G entropy + spectral β + perturbation weight into a robust 3D Class IV discriminator?
- The 770 "tracks" in Rule 30 — can we develop a false-positive rejection criterion?
- Connection to Brooklyn's SAT phase transition: is the commutator measuring constraint density?

## Heartbeat Protocol
Cycles alternate: **BUILD** (add experiments, run code) → **AUDIT** (review claims, check rigor, prune speculation).

### BUILD cycle tasks (pick one per cycle):
1. Try improved spectral methods (Welch, DFA) for cleaner β estimates
2. Run larger scaling experiments (N=1001+) for Rule 54
3. Investigate Class II high-entropy outliers (Rules 24, 35)
4. **D-as-input experiment**: Feed D(S) back as input to the rule — CA with self-awareness of its own dynamics. Does this produce fixed points, novel oscillations, or emergent memory?
5. **Novel operator compositions**: Explore D(I(S,d)), C(D(S), D₂(S)), I(S, G(S)), etc. Look for replicable structure that "shouldn't make sense."
6. **Multi-feature classifier**: Combine G entropy + spectral β + perturbation weight + nested ratio for Class IV prediction
7. Extend to 2D CA or larger rule spaces

### AUDIT cycle tasks:
1. Review all conjectures — are they supported by the corrected data?
2. Check for overclaiming in the paper
3. Verify experimental methodology (sample sizes, statistical significance)
4. Ensure PLAN.md is current
5. Clean up code quality, add docstrings, ensure reproducibility

### North Star
**An algorithmic Class IV detector.** Everything serves this goal.

### Guiding Principle
D, I, C all return configuration spaces of the same dimensionality as the input. This means arbitrary composition is well-typed. Explore the operator algebra freely — novel compositions and playful exploration may be rewarded. We're looking for replicable behavior that shouldn't make sense.

### Always:
- Update this PLAN.md when done
- Git commit and push
- Report significant findings to Discord #general
