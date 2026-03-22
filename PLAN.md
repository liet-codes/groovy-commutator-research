# PLAN.md — Current Status and Next Steps

*Last updated: 2026-03-22*

## Completed
- [x] Core library: ECA engine, operators, measures, affine decomposition
- [x] Experiment 1: Baseline characterization (all 256 rules)
- [x] Experiment 2: Commutator zoo (8 variants × 20 rules)
- [x] Experiment 3: Class IV hunt (extended runs)
- [x] Experiment 4: Affine × commutator analysis
- [x] Paper draft v1: definitions, results, conjectures
- [ ] Experiment 5: 88 equivalence class lookup table (WRITTEN, NOT YET RUN)

## Key Findings So Far
1. Class IV rules have highest mean G entropy (0.892) with tightest std (0.079)
2. Class IV has longest correlation lengths (1.63 mean)
3. No Class IV rule is affine (0% vs 27% of Class III)
4. Nested commutator [E,[E,D]] drops significantly for Rule 110 (0.973 → 0.713)
5. Multi-step [Eⁿ,D] grows for Class IV but saturates for Class III
6. Conjecture: Class IV = non-affine AND G entropy > 0.8 AND G corr length > 1.3

## Immediate Next Steps
1. **Run exp5_equivalence_table.py** — generates the 88-class lookup table (Brooklyn needs this)
2. **Push to GitHub** as `liet-codes/groovy-commutator-research` (public)
3. **Add equivalence table to paper** as Appendix A

## Medium-Term
4. **Spectral analysis of G time series** — compute power spectra, look for 1/f signatures
5. **Glider detection via G** — do local G peaks correspond to glider positions?
6. **Larger systems** — run key experiments at N=501, N=1001 to check finite-size effects
7. **2D CA** — extend commutator analysis to 2-state 2D totalistic rules

## Open Questions
- Why does the nested commutator drop for Rule 110 but not Rule 54? Is there sub-structure within Class IV?
- Is there a commutator variant that perfectly separates Class IV? The zoo results suggest [E,[E,D]] might be key.
- Can we prove Conjecture 1 (commutator vanishing for affine rules) algebraically?
- Does the 1/9 universal tiling ratio (Patrick's result) show up in commutator statistics?
- Connection to Brooklyn's SAT phase transition: is the commutator measuring constraint density?

## Heartbeat Tasks
When running on heartbeat cycles:
1. Check if any experiments need re-running with updated parameters
2. Run exp5 if not yet complete
3. Try new commutator variants suggested in RESEARCH.md
4. Update paper with any new results
5. Always update this PLAN.md when done
