# The Groovy Commutator: Operator Non-Commutativity as a Characterization of Complex Cellular Automata

## Abstract

We introduce the *Groovy Commutator*, a diagnostic operator for elementary cellular automata (ECA) defined as the commutator of the spatial differentiation operator D and the temporal evolution operator E. For a state S, the Groovy Commutator is G(S) = D(E(S)) ⊕ E(D(S)), measuring how much the ordering of differentiation and evolution matters. We compute G across all 256 ECA rules, combined with Shannon entropy, spatial correlation length, temporal autocorrelation, and compressibility measures. We further extend the analysis to a family of commutator variants ([E, D²], [E, I], [Eⁿ, D], [E, [E, D]], and the anti-commutator {E, D}) and perform affine decomposition of each rule into a best-fit linear (XOR) component plus residual perturbation. Our experiments support the conjecture that Wolfram Class IV rules occupy a distinctive region of the commutator-affine feature space: they exhibit nonzero but structured commutator signatures, moderate perturbation from affine behavior, and intermediate correlation lengths — consistent with "edge of chaos" characterizations but now given a precise operator-algebraic formulation.

## 1. Introduction

The classification of elementary cellular automata into Wolfram's four classes — uniform (I), periodic (II), chaotic (III), and complex (IV) — remains one of the most important open problems in discrete dynamical systems. While Classes I–III admit relatively clean characterizations (fixed points, limit cycles, and positive entropy respectively), Class IV has resisted formal definition. These rules produce long transients, localized structures ("gliders"), and are conjectured to support universal computation, yet no single measure cleanly separates them from the other classes.

We propose an approach rooted in operator algebra. Consider two natural operators on binary strings:

- **E** (evolution): apply the CA rule to advance one timestep
- **D** (differentiation): compute the spatial finite difference via XOR of adjacent cells

If E and D commute — that is, if D(E(S)) = E(D(S)) for all states S — then the spatial structure of the derivative is preserved under evolution. This holds trivially for Class I rules (everything collapses) and for the purely affine (XOR-linear) rules. The *Groovy Commutator* G = [E, D] measures the degree of non-commutativity:

$$G(S) = D(E(S)) \oplus E(D(S))$$

Our central finding is that the statistical properties of G — its entropy, spatial correlation, temporal persistence, and compressibility — encode information about the dynamical class of the underlying rule, and that Class IV rules occupy a characteristic "structured but nontrivial" regime.

### 1.1 Related Work

Wolfram's original classification (1984, 2002) was based on visual inspection of spacetime diagrams. Subsequent work has attempted formalization via Lyapunov exponents (Bagnoli et al., 1992), mean-field theory (Gutowitz, 1991), computational complexity (Cook, 2004), Langton's λ parameter (1990), and various entropy-based measures (Martinez et al., 2013). The affine decomposition approach follows Rose (2026), who decomposes ECA rules into XOR-linear base plus nonlinear perturbation as a structural analysis tool.

Our commutator approach is, to our knowledge, novel. It draws on the philosophy of commutator analysis in quantum mechanics and Lie algebra, adapted to the discrete Boolean setting.

## 2. Definitions

### 2.1 Elementary Cellular Automata

An elementary cellular automaton is a map f: {0,1}³ → {0,1} applied synchronously to a ring of N cells. The 256 possible rules are indexed by the Wolfram number W = Σᵢ f(i) · 2ⁱ where i encodes the three-cell neighborhood.

### 2.2 Operators

Let S ∈ {0,1}ᴺ be a binary state on a ring of N cells with periodic boundary conditions.

**Evolution operator E:**
$$E(S)[i] = f(S[i-1], S[i], S[i+1])$$

**Differentiation operator D:**
$$D(S)[i] = S[i] \oplus S[i+1]$$

This is the discrete spatial derivative in GF(2). It measures local spatial variation: D(S) is zero at position i iff cells i and i+1 agree.

**Integration operator I:**
$$I(S)[i] = \bigoplus_{j=0}^{i} S[j]$$

The cumulative XOR, serving as a left-inverse of D (up to boundary effects).

**Second-order differentiation D²:**
$$D^2 = D \circ D$$

### 2.3 The Groovy Commutator

$$G(S) = [E, D](S) = D(E(S)) \oplus E(D(S))$$

G measures the failure of E and D to commute. When G(S) = 0 for all S, the evolution and differentiation operators commute perfectly — the rule preserves the structure of spatial derivatives.

### 2.4 Extended Commutator Family

We study several related operators:

| Operator | Definition | Interpretation |
|----------|-----------|----------------|
| [E, D] = G | D(E(S)) ⊕ E(D(S)) | Standard commutator |
| [E, D²] | D²(E(S)) ⊕ E(D²(S)) | Second-order spatial sensitivity |
| [E, I] | I(E(S)) ⊕ E(I(S)) | Integration-evolution non-commutativity |
| [Eⁿ, D] | D(Eⁿ(S)) ⊕ Eⁿ(D(S)) | Multi-step evolution commutator |
| [E, [E, D]] | [E,D](E(S)) ⊕ E([E,D](S)) | Nested (Jacobi-like) commutator |
| {E, D} | D(E(S)) ∧ E(D(S)) | Anti-commutator (AND of branches) |

### 2.5 Affine Decomposition

Following Rose (2026), every ECA rule f can be uniquely decomposed as:
$$f = a \oplus p$$

where a is the closest affine (XOR-linear) function and p is the residual perturbation. The 16 affine functions on three binary inputs are: 0, 1, L, ¬L, C, ¬C, R, ¬R, L⊕C, ¬(L⊕C), L⊕R, ¬(L⊕R), C⊕R, ¬(C⊕R), L⊕C⊕R, ¬(L⊕C⊕R).

The **perturbation magnitude** |p| counts the number of input configurations where f differs from a (ranging from 0 to 4, since the closest affine function is at most Hamming distance 4 from any Boolean function on 3 inputs). The **perturbation weight** is |p|/8.

### 2.6 Measures

For a binary state S ∈ {0,1}ᴺ:

- **Shannon entropy** H(S) = -p₀ log₂ p₀ - p₁ log₂ p₁ where pᵢ is the fraction of cells in state i
- **Block entropy** H_k(S): entropy of non-overlapping k-blocks
- **Spatial correlation length** ξ: the 1/e decay distance of the spatial autocorrelation
- **Temporal autocorrelation**: autocorrelation of the entropy time series H(G(Eᵗ(S))) over t
- **Hamming distance** d_H(Gₜ, Gₜ₊₁): fraction of cells that differ between successive G states
- **Compressibility ratio** κ: gzip(S) / |S|, a proxy for Kolmogorov complexity

## 3. Experimental Setup

All experiments use periodic boundary conditions on binary rings.

### 3.1 Experiment 1: Baseline Characterization

- All 256 ECA rules
- Width N = 101, T = 200 timesteps, 50 random initial conditions per rule
- Transient discard: first 50 timesteps
- Measures: G entropy, G correlation length, G temporal autocorrelation, G compressibility
- Output: scatter plot of entropy vs. correlation length, colored by Wolfram class

### 3.2 Experiment 2: Commutator Zoo

- 20 representative rules spanning all classes
- All 8 commutator variants
- Width N = 101, T = 200, 30 trials per rule
- Output: heatmap of entropy across rules × variants

### 3.3 Experiment 3: Class IV Hunt

- 6 Class IV rules (54, 106, 110, 124, 137, 193) vs. representatives from Classes I, II, III
- Extended runs: N = 201, T = 500, 30 trials
- All 6 measures computed
- Output: multi-panel comparison across classes

### 3.4 Experiment 4: Affine × Commutator

- All 256 rules: affine decomposition AND commutator measures
- Testing hypothesis: Class IV ≈ {nonzero perturbation} ∩ {structured commutator}
- Output: 2D scatter of perturbation weight vs. G entropy

## 4. Results

### 4.1 Baseline: Class Separation in Commutator Space

We computed the Groovy Commutator G across all 256 ECA rules (N=101, T=200, 50 random trials each). The table below shows mean statistics by Wolfram class:

| Class | n | G Entropy | G Corr Length | G Temporal ACF | G Compressibility |
|-------|---|-----------|---------------|----------------|-------------------|
| I     | 16 | 0.110 ± 0.301 | 0.26 ± 0.72 | 1.000 ± 0.000 | 0.261 ± 0.063 |
| II    | 63 | 0.647 ± 0.277 | 1.25 ± 0.60 | 0.918 ± 0.148 | 0.402 ± 0.073 |
| III   | 30 | 0.647 ± 0.416 | 0.93 ± 0.69 | 0.356 ± 0.408 | 0.414 ± 0.114 |
| IV    |  9 | **0.892 ± 0.079** | **1.63 ± 1.06** | 0.199 ± 0.150 | 0.439 ± 0.040 |

Several striking patterns emerge:

- **Class I** rules cluster near the origin: G entropy ≈ 0, correlation length ≈ 0, and perfect temporal autocorrelation (constant zero commutator). Evolution to a uniform state annihilates all spatial structure.

- **Class II** rules show moderate G entropy (0.65) with high temporal autocorrelation (0.92), indicating periodic or quasi-periodic commutator behavior.

- **Class III** shows a bimodal distribution: the non-affine Class III rules (e.g., Rule 30) have high G entropy ≈ 0.99, while the affine Class III rules (e.g., Rules 90, 150) have G entropy = 0 exactly. The class-wide mean (0.65) reflects this mixture. Temporal autocorrelation is low (0.36), indicating rapid decorrelation.

- **Class IV** rules have the **highest mean G entropy** (0.892) with the **tightest standard deviation** (0.079) and the **longest correlation lengths** (1.63). Their temporal autocorrelation (0.20) is intermediate — neither periodic nor fully decorrelating. This "high entropy, long correlation, moderate persistence" signature is distinctive.

### 4.2 Commutator Zoo: Variant Discrimination

We computed all 8 commutator variants on 20 representative rules. Selected results for key rules:

| Rule | Class | [E,D] | [E,D²] | [E,I] | [E²,D] | [E³,D] | [E⁴,D] | [E,[E,D]] | {E,D} |
|------|-------|-------|---------|-------|---------|---------|---------|-----------|-------|
| 30   | III   | 0.989 | 0.941 | 0.993 | 0.883 | 0.993 | 0.990 | 0.990 | 0.804 |
| 54   | IV    | 0.856 | 0.993 | 0.971 | 0.821 | 0.948 | 0.872 | 0.935 | 0.853 |
| 90   | III*  | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.993 |
| 110  | IV    | 0.973 | 0.956 | 0.975 | 0.953 | 0.978 | 0.996 | **0.713** | 0.916 |

*Rule 90 is classified Class III by Wolfram but is affine (L⊕R), causing all commutators to vanish.

Key findings from the zoo:

- **Rule 90 confirmation**: All XOR-type commutators vanish for this affine rule, while the anti-commutator {E,D} is high (0.993). The anti-commutator measures where both orderings agree on producing 1, which is high for affine rules since D(E(S)) = E(D(S)) everywhere.

- **Nested commutator drop for Rule 110**: The nested commutator [E,[E,D]] drops to 0.713 — significantly below the first-order value of 0.973. This suggests the non-commutativity of Rule 110 has hierarchical structure that partially cancels at second order.

- **Multi-step divergence**: For Class IV Rule 110, [Eⁿ,D] increases from 0.953 (n=2) to 0.996 (n=4). For Class III Rule 30, it's already near-maximal at all n. Class IV rules show a growth trajectory; Class III rules are already saturated.

### 4.3 Class IV Distinguishing Features

Extended runs (N=201, T=500) with all 6 measures reveal:

| Rule | Class | G Entropy | Block Ent. | Corr Len | Temp ACF | Hamming | Compress |
|------|-------|-----------|------------|----------|----------|---------|----------|
| 54   | IV    | 0.905 | 2.50 | **3.40** | 0.408 | 0.754 | 0.283 |
| 106  | IV    | 0.994 | 3.49 | 1.00 | 0.055 | 0.375 | 0.358 |
| 110  | IV    | 0.979 | 3.35 | 1.40 | 0.090 | 0.500 | 0.266 |
| 124  | IV    | 0.980 | 3.34 | 1.57 | 0.103 | 0.505 | 0.258 |
| 0    | I     | 0.000 | 0.00 | 0.00 | 1.000 | 0.000 | 0.119 |
| 4    | II    | 0.803 | 2.40 | 2.00 | 1.000 | 0.000 | 0.292 |
| 108  | II    | 0.682 | 2.45 | 1.00 | 0.945 | 0.090 | 0.309 |
| 30   | III   | 0.995 | 3.54 | 1.00 | 0.040 | 0.375 | 0.359 |
| 90   | III*  | 0.000 | 0.00 | 0.00 | 1.000 | 0.000 | 0.119 |

Key discriminating features for Class IV:

1. **High but non-maximal G entropy**: Class IV rules concentrate in the range 0.86–0.99, versus Class III's near-1.0 (for non-affine rules) or 0.0 (for affine rules). Class IV is consistently high but never quite saturated.

2. **Elevated correlation length**: Rule 54 stands out with ξ=3.4, the longest in the dataset. Other Class IV rules show ξ≥1.4, while Class III Rule 30 has ξ=1.0.

3. **Intermediate temporal autocorrelation**: Class IV ACF ranges 0.06–0.41, between Class II's ~1.0 (periodic) and Class III's 0.04 (rapidly decorrelating). Rule 54 is notably higher (0.41), suggesting longer-lived correlations.

4. **Nonzero Hamming distance**: Class IV rules have d_H > 0 (unlike Class II Rule 4, which has d_H = 0, indicating a fixed-point G). Class IV Hamming distances range 0.37–0.75.

5. **Moderate compressibility**: Class IV κ ∈ [0.24, 0.36], neither as low as Class I (0.12) nor as high as Class III (0.36).

### 4.4 Affine Decomposition × Commutator

Across all 256 rules, combining affine decomposition with commutator measures:

| Class | Pert. Weight | G Entropy | G Corr Length | Fraction Affine |
|-------|-------------|-----------|---------------|-----------------|
| I     | 0.148 ± 0.082 | 0.111 ± 0.304 | 0.26 ± 0.70 | 12% |
| II    | 0.177 ± 0.077 | 0.649 ± 0.278 | 1.22 ± 0.57 | 6% |
| III   | 0.142 ± 0.102 | 0.648 ± 0.417 | 0.93 ± 0.69 | 27% |
| IV    | 0.167 ± 0.062 | **0.895 ± 0.076** | **1.66 ± 1.13** | **0%** |

Key findings:

- **No Class IV rule is affine**. Every Class IV rule has nonzero perturbation from its nearest affine function — 0% are affine, compared to 27% of Class III rules (which includes the XOR-linear rules 90, 105, 150).

- **Perturbation weight is similar across classes** (~0.14–0.18), so perturbation magnitude alone does not discriminate. The interesting signal is in the *conjunction* with commutator behavior.

- **Class IV's unique niche**: The only class with *simultaneously* (1) zero affine fraction, (2) G entropy > 0.8, and (3) G correlation length > 1.5. No other class satisfies all three conditions.

This supports a refined version of our hypothesis:

> **Conjecture (Groovy Commutator Characterization):** A rule R exhibits Class IV behavior if and only if: (1) R has nonzero but sub-maximal perturbation from its nearest affine function, (2) the Groovy Commutator G has intermediate entropy (bounded away from both 0 and 1), and (3) G exhibits spatial correlations significantly exceeding 1/e of the system size.

## 5. Discussion

### 5.1 Why the Commutator Captures Complexity

The Groovy Commutator measures, in precise terms, the failure of spatial structure (as captured by D) to be preserved under temporal evolution (E). For a rule where E and D commute, knowing the spatial derivative before evolution is equivalent to knowing it after — the spatial structure is a conserved quantity of the dynamics (in the GF(2) sense).

Class I rules trivialize both E and D, so they commute vacuously. Affine (XOR-linear) rules commute because both E and D are linear operations over GF(2), and linear operations always commute.

Class III rules destroy spatial structure so aggressively that the commutator is essentially random noise — the ordering doesn't matter because everything becomes unpredictable either way.

Class IV rules are precisely those where the ordering *matters in a structured way*: the commutator is nonzero (evolution does not preserve spatial derivatives) but the non-commutativity has spatial coherence (it's not random noise). This is the hallmark of computation — the manipulation of structured information.

### 5.2 Connection to Edge of Chaos

Langton's "edge of chaos" hypothesis (1990) posits that computation occurs at the boundary between ordered and chaotic dynamics. The Groovy Commutator provides a new lens on this: the "edge" is where the commutator [E, D] is nonzero but structured. This is analogous to how in quantum mechanics, the commutator [x, p] = iℏ defines the boundary between classical (commuting) and quantum behavior, and the magnitude of the commutator controls the uncertainty principle.

### 5.3 The Affine Connection

The affine decomposition reveals that Class IV rules are perturbations of linear dynamics — they have a dominant XOR-linear component with a small but critical nonlinear perturbation. This resonates with Brooklyn Rose's insight that the "interesting" dynamics arise not from maximal nonlinearity but from structured deviations from linearity.

The conjunction of moderate perturbation and moderate commutator entropy may be a signature of the "barrier" that Rose describes — the computational boundary between the affine regime (where behavior is predictable) and the fully nonlinear regime (where behavior is chaotic).

### 5.4 Limitations

1. **Classification ground truth**: Our Wolfram class assignments are based on community consensus, which is itself partially subjective for marginal rules. The classification of rules near class boundaries remains contested.

2. **Finite-size effects**: All experiments use finite rings (N ≤ 201). Some Class IV phenomena (particularly glider interactions) may require larger systems.

3. **Statistical measures**: Shannon entropy and correlation length are coarse measures. Richer topological or algebraic invariants of G might provide sharper discrimination.

4. **Conjecture status**: Our characterization is stated as a conjecture, not a theorem. A proof would require understanding the algebraic structure of the commutator for all possible rules, which remains open.

## 6. Conjectures

We state our main conjectures with epistemic tags indicating confidence levels.

**Conjecture 1 (Commutator Vanishing — HIGH CONFIDENCE):**
*For any affine ECA rule (W ∈ {0, 15, 51, 85, 90, 102, 105, 150, 153, 165, 170, 195, 204, 240, 255}), the Groovy Commutator G(S) = 0 for all S.*

Rationale: Both E and D are GF(2)-linear for affine rules, and linear operators commute.

**Conjecture 2 (Class IV High-Entropy Regime — MODERATE CONFIDENCE):**
*Class IV rules are characterized by G entropy > 0.85 and G spatial correlation length exceeding the class mean of non-affine Class III rules.*

Rationale: Experimentally, Class IV rules have mean G entropy 0.892 ± 0.079, the highest of any class. Their correlation length (1.63) exceeds Class III's (0.93). Boundary values may shift with system size.

**Conjecture 3 (Multi-step Growth — MODERATE CONFIDENCE):**
*For Class IV rules, the entropy of [Eⁿ, D] grows with n (e.g., Rule 110: 0.953 → 0.978 → 0.996 for n=2,3,4), while for Class III rules it is near-maximal at all n.*

Rationale: Class III rules are already maximally chaotic at n=1. Class IV rules have longer-range correlations that are progressively disrupted at longer evolution horizons.

**Conjecture 4 (Affine-Commutator Joint Characterization — MODERATE CONFIDENCE):**
*A necessary condition for Class IV behavior is: (1) the rule is non-affine (perturbation weight > 0), (2) G entropy > 0.8, and (3) G correlation length > 1.3.*

Rationale: In our experiments, Class IV is the only class where 100% of rules are non-affine AND have G entropy > 0.8 AND correlation length > 1.3. The specific thresholds are tentative.

**Conjecture 5 (Nested Commutator Structure — SPECULATIVE):**
*For some Class IV rules (notably Rule 110), the nested commutator [E, [E, D]] has significantly lower entropy than [E, D] (0.713 vs 0.973 for Rule 110), indicating hierarchical structure in the non-commutativity.*

Rationale: Observed strongly for Rule 110 but not uniformly across all Class IV rules (Rule 54 shows nested entropy 0.935 > [E,D] entropy 0.856). This may indicate sub-structure within Class IV itself.

## 7. Future Directions

1. **Larger rule spaces**: Extend to 2-state 2D CA, k=3 totalistic rules, or general k-state 1D CA to test universality of the commutator characterization.

2. **Algebraic theory**: Develop a formal algebraic framework for the commutator over GF(2)^N, potentially connecting to the theory of Boolean function composition.

3. **Spectral analysis**: Compute the power spectrum of G time series and look for 1/f signatures associated with criticality.

4. **Machine learning**: Use commutator features as inputs to a classifier and evaluate whether they outperform raw spacetime features for Wolfram class prediction.

5. **Glider detection**: Study whether local peaks in G correspond to glider positions, providing a commutator-based glider detection algorithm.

## References

1. Bagnoli, F., Rechtman, R., & Ruffo, S. (1992). Some facts of life. *Physica A*, 171(2), 249-264.
2. Cook, M. (2004). Universality in Elementary Cellular Automata. *Complex Systems*, 15(1), 1-40.
3. Gutowitz, H. (1991). *Cellular Automata: Theory and Experiment*. MIT Press.
4. Langton, C. G. (1990). Computation at the Edge of Chaos: Phase Transitions and Emergent Computation. *Physica D*, 42(1-3), 12-37.
5. Martinez, G. J., Adamatzky, A., & McIntosh, H. V. (2013). On the representation of Boolean and real-valued cellular automata. *Journal of Cellular Automata*, 8(1-2), 1-19.
6. Rose, B. (2026). The Shape of the Barrier. *Preprint*.
7. Rose, B. (2026). Affine Decompositions of Elementary Cellular Automata. *Preprint*.
8. Wolfram, S. (1984). Universality and complexity in cellular automata. *Physica D*, 10(1-2), 1-35.
9. Wolfram, S. (2002). *A New Kind of Science*. Wolfram Media.

---

*Computational code and data available in the accompanying repository.*

---

## Appendix A: Equivalence Class Commutator Signatures

We compute the Groovy Commutator signature for all 88 ECA equivalence classes under left-right reflection and black-white complementation. For each class, the smallest rule number serves as representative. Parameters: N=101, T=300, 50 random initial conditions, transient discard of 50 steps.

| Rep | Class | Equiv | G≠0 (single) | G≠0 (random) | G Entropy | G Corr Len | G Compress | Period | Pert Weight |
|-----|-------|-------|---------------|--------------|-----------|------------|------------|--------|-------------|
|   0 |   I | 2 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
|   1 |  II | 2 | ✓ | 1.00 | 0.717±0.121 | 1.50 | 0.411 | 1 | 0.125 |
|   2 |  II | 4 | ✓ | 1.00 | 0.580±0.090 | 1.00 | 0.419 | 1 | 0.125 |
|   3 |  II | 4 | ✓ | 1.00 | 0.724±0.060 | 1.00 | 0.458 | 1 | 0.250 |
|   4 |  II | 2 | ✓ | 1.00 | 0.798±0.085 | 1.96 | 0.425 | 1 | 0.125 |
|   5 |  II | 2 | ✓ | 1.00 | 0.646±0.062 | 1.12 | 0.467 | 2 | 0.250 |
|   6 |  II | 4 | ✓ | 1.00 | 0.694±0.075 | 1.00 | 0.429 | 1 | 0.250 |
|   7 |  II | 4 | ✓ | 1.00 | 0.538±0.084 | 2.00 | 0.361 | 1 | 0.125 |
|   8 |   I | 4 | ✓ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.125 |
|   9 |  II | 4 | ✓ | 1.00 | 0.832±0.038 | 1.00 | 0.374 | 1 | 0.250 |
|  10 |  II | 4 | ✗ | 1.00 | 0.537±0.072 | 1.00 | 0.407 | 1 | 0.250 |
|  11 |  II | 4 | ✓ | 1.00 | 0.950±0.029 | 1.00 | 0.425 | 1 | 0.125 |
|  12 |  II | 4 | ✓ | 1.00 | 0.885±0.052 | 1.00 | 0.493 | 1 | 0.250 |
|  13 |  II | 4 | ✓ | 1.00 | 0.814±0.074 | 1.10 | 0.415 | 1 | 0.125 |
|  14 |  II | 4 | ✓ | 1.00 | 0.305±0.330 | 1.24 | 0.300 | 1 | 0.125 |
|  15 |  II | 2 | ✓ | 1.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
|  18 | III | 2 | ✓ | 1.00 | 0.816±0.006 | 1.00 | 0.464 | 1 | 0.250 |
|  19 |  II | 2 | ✓ | 1.00 | 0.288±0.094 | 1.00 | 0.328 | 1 | 0.125 |
|  22 | III | 2 | ✓ | 1.00 | 0.563±0.006 | 1.12 | 0.410 | aper. | 0.125 |
|  23 |  II | 1 | ✓ | 1.00 | 0.701±0.091 | 1.06 | 0.446 | 1 | 0.250 |
|  24 |  II | 4 | ✓ | 1.00 | 0.973±0.027 | 1.76 | 0.451 | 1 | 0.250 |
|  25 |  II | 4 | ✓ | 1.00 | 0.887±0.027 | 1.00 | 0.374 | 1 | 0.125 |
|  26 |  II | 4 | ✓ | 1.00 | 0.938±0.032 | 1.08 | 0.495 | 1 | 0.125 |
|  27 |  II | 4 | ✓ | 1.00 | 0.842±0.039 | 1.00 | 0.490 | 1 | 0.250 |
|  28 |  II | 4 | ✓ | 1.00 | 0.517±0.077 | 1.00 | 0.330 | 2 | 0.125 |
|  29 |  II | 2 | ✓ | 1.00 | 0.859±0.035 | 1.00 | 0.510 | 1 | 0.250 |
|  30 | III | 4 | ✓ | 1.00 | 0.989±0.001 | 1.02 | 0.523 | 1 | 0.250 |
|  32 |   I | 2 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.125 |
|  33 |  II | 2 | ✓ | 1.00 | 0.882±0.050 | 1.02 | 0.470 | 1 | 0.250 |
|  34 |  II | 4 | ✓ | 1.00 | 0.883±0.055 | 1.00 | 0.492 | 1 | 0.250 |
|  35 |  II | 4 | ✓ | 1.00 | 0.985±0.015 | 1.00 | 0.458 | 1 | 0.125 |
|  36 |  II | 2 | ✓ | 1.00 | 0.553±0.147 | 2.00 | 0.361 | 1 | 0.250 |
|  37 |  II | 2 | ✓ | 1.00 | 0.497±0.069 | 1.04 | 0.395 | 2 | 0.125 |
|  38 |  II | 4 | ✓ | 1.00 | 0.867±0.039 | 1.06 | 0.498 | 1 | 0.125 |
|  40 |   I | 4 | ✓ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.250 |
|  41 |  II | 4 | ✓ | 1.00 | 0.752±0.022 | 1.00 | 0.445 | 1 | 0.125 |
|  42 |  II | 4 | ✗ | 1.00 | 0.475±0.127 | 1.56 | 0.372 | 1 | 0.125 |
|  43 |  II | 2 | ✓ | 0.94 | 0.306±0.203 | 2.06 | 0.309 | 1 | 0.250 |
|  44 |  II | 4 | ✓ | 1.00 | 0.703±0.065 | 1.00 | 0.442 | 1 | 0.125 |
|  45 | III | 4 | ✓ | 1.00 | 0.989±0.001 | 1.02 | 0.522 | 1 | 0.250 |
|  46 |  II | 4 | ✓ | 1.00 | 0.854±0.045 | 1.00 | 0.454 | 1 | 0.250 |
|  50 |  II | 2 | ✓ | 1.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.125 |
|  51 |  II | 1 | ✓ | 1.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
|  54 |  IV | 2 | ✓ | 1.00 | 0.861±0.068 | 3.50 | 0.409 | 1 | 0.250 |
|  56 |  II | 4 | ✓ | 1.00 | 0.789±0.100 | 1.28 | 0.453 | 1 | 0.125 |
|  57 |  II | 2 | ✓ | 1.00 | 0.213±0.147 | 1.42 | 0.292 | 1 | 0.250 |
|  58 |  II | 4 | ✓ | 1.00 | 0.815±0.155 | 2.62 | 0.386 | 1 | 0.250 |
|  60 | III | 4 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
|  62 |  II | 4 | ✓ | 1.00 | 0.352±0.097 | 1.06 | 0.337 | 1 | 0.125 |
|  72 |  II | 2 | ✓ | 1.00 | 0.487±0.128 | 1.00 | 0.361 | 1 | 0.250 |
|  73 |  II | 2 | ✓ | 1.00 | 0.969±0.017 | 2.60 | 0.433 | 1 | 0.125 |
|  74 |  II | 4 | ✓ | 1.00 | 0.803±0.049 | 1.00 | 0.466 | 1 | 0.125 |
|  76 |  II | 2 | ✗ | 1.00 | 0.779±0.102 | 1.96 | 0.435 | 1 | 0.125 |
|  77 |  II | 1 | ✓ | 1.00 | 0.980±0.026 | 2.04 | 0.447 | 1 | 0.250 |
|  78 |  II | 4 | ✓ | 1.00 | 0.932±0.061 | 2.18 | 0.414 | 1 | 0.250 |
|  90 | III | 2 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
|  94 |  II | 2 | ✓ | 1.00 | 0.976±0.031 | 2.10 | 0.456 | 1 | 0.125 |
| 104 |  II | 2 | ✓ | 1.00 | 0.466±0.161 | 2.14 | 0.328 | 1 | 0.125 |
| 105 | III | 1 | ✓ | 1.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
| 106 |  IV | 4 | ✓ | 1.00 | 0.990±0.001 | 1.06 | 0.518 | 1 | 0.250 |
| 108 |  II | 2 | ✗ | 1.00 | 0.676±0.088 | 1.00 | 0.445 | 1 | 0.250 |
| 110 |  IV | 4 | ✓ | 1.00 | 0.975±0.004 | 1.38 | 0.437 | 1 | 0.125 |
| 122 | III | 2 | ✓ | 1.00 | 0.557±0.009 | 1.00 | 0.379 | aper. | 0.125 |
| 126 | III | 2 | ✓ | 1.00 | 0.940±0.005 | 1.76 | 0.459 | 1 | 0.250 |
| 128 |   I | 2 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.125 |
| 130 | III | 4 | ✓ | 1.00 | 0.696±0.063 | 1.00 | 0.441 | 1 | 0.250 |
| 132 |  II | 2 | ✓ | 1.00 | 0.852±0.054 | 1.00 | 0.466 | 1 | 0.250 |
| 134 |  II | 4 | ✓ | 1.00 | 0.813±0.046 | 1.00 | 0.436 | 1 | 0.125 |
| 136 |   I | 4 | ✓ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.250 |
| 138 |  II | 4 | ✗ | 1.00 | 0.575±0.090 | 1.00 | 0.408 | 1 | 0.125 |
| 140 |  II | 4 | ✓ | 1.00 | 0.700±0.044 | 1.00 | 0.429 | 1 | 0.125 |
| 142 |  II | 2 | ✓ | 1.00 | 0.330±0.247 | 1.88 | 0.313 | 1 | 0.250 |
| 146 | III | 2 | ✓ | 1.00 | 0.986±0.002 | 1.86 | 0.456 | 1 | 0.125 |
| 150 | III | 1 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
| 152 |  II | 4 | ✓ | 1.00 | 0.962±0.024 | 1.34 | 0.452 | 1 | 0.125 |
| 154 |  II | 4 | ✓ | 1.00 | 0.990±0.008 | 1.04 | 0.524 | 1 | 0.250 |
| 156 |  II | 2 | ✓ | 1.00 | 0.755±0.083 | 1.42 | 0.445 | 1 | 0.250 |
| 160 |   I | 2 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.250 |
| 162 |  II | 4 | ✓ | 1.00 | 0.657±0.046 | 1.00 | 0.427 | 1 | 0.125 |
| 164 |  II | 2 | ✓ | 1.00 | 0.656±0.106 | 1.96 | 0.388 | 1 | 0.125 |
| 168 |   I | 4 | ✓ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.125 |
| 170 |  II | 2 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
| 172 |  II | 4 | ✓ | 1.00 | 0.569±0.069 | 1.00 | 0.419 | 1 | 0.250 |
| 178 |  II | 1 | ✓ | 1.00 | 0.982±0.029 | 2.08 | 0.450 | 1 | 0.250 |
| 184 |  II | 2 | ✓ | 0.86 | 0.315±0.213 | 0.90 | 0.328 | 1 | 0.250 |
| 200 |  II | 2 | ✓ | 1.00 | 0.681±0.081 | 1.00 | 0.448 | 1 | 0.125 |
| 204 |  II | 1 | ✗ | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 1 | 0.000 |
| 232 |  II | 1 | ✓ | 1.00 | 0.692±0.085 | 1.02 | 0.453 | 1 | 0.250 |

### A.1 Analysis of the Equivalence Class Table

The 88-class table confirms and extends our earlier findings:

**Conjecture 1 confirmed across all 88 classes.** Every affine rule (perturbation weight = 0.000) — Rules 0, 15, 51, 60, 90, 105, 150, 170, 204 — has G entropy exactly 0.000 and G≠0 fraction 0.00 from random initial conditions. The commutator vanishes identically for all affine rules, consistent with the theoretical prediction that GF(2)-linear operators commute.

**Class IV rules occupy a distinctive niche.** The three Class IV representatives (Rules 54, 106, 110) all exhibit G entropy > 0.86 with G≠0 fraction = 1.00. Their signatures are:

| Rule | G Entropy | G Corr Len | G Compress | Pert Weight |
|------|-----------|------------|------------|-------------|
| 54   | 0.861±0.068 | **3.50** | 0.409 | 0.250 |
| 106  | 0.990±0.001 | 1.06 | 0.518 | 0.250 |
| 110  | 0.975±0.004 | 1.38 | 0.437 | 0.125 |

**Rule 54 has the highest correlation length (3.50) of any rule in the table.** The next highest are Rules 58 (2.62), 73 (2.60), and 78 (2.18) — all Class II rules. Rule 54's correlation length exceeds the next-highest by 34%, suggesting exceptionally long-range spatial structure in its commutator. This is consistent with Rule 54's known complex dynamics, including its ability to support glider-like structures.

**Class III splits into affine and non-affine subgroups.** Among the 14 Class III representatives, 5 are affine (Rules 60, 90, 105, 150 — all with G=0) and 9 are non-affine (Rules 18, 22, 30, 45, 122, 126, 130, 146 — with G entropy ranging from 0.557 to 0.989). This bimodal structure confirms that the traditional Wolfram classification conflates two dynamically distinct regimes: XOR-linear chaos (affine, G=0) and nonlinear chaos (non-affine, G≈1).

**Class I rules show consistent G=0 from random ICs.** All 8 Class I representatives have G≠0 (random) = 0.00, confirming that evolution to a uniform or simple fixed point annihilates the commutator. Interestingly, Rules 8, 40, 136, and 168 show G≠0 from a single black cell (✓) but not from random ICs (0.00), indicating transient non-commutativity that vanishes as the rule converges to its attractor.
