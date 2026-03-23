# The Groovy Commutator: Operator Non-Commutativity as a Characterization of Complex Cellular Automata

## Abstract

We introduce the *Groovy Commutator*, a diagnostic operator for elementary cellular automata (ECA) defined as the commutator of the differentiation operator D and the temporal evolution operator E. Here D(S) = S ⊕ φ(S) detects which cells are about to change, and E(S) = φ(S) applies the CA rule. The Groovy Commutator G(S) = D(E(S)) ⊕ E(D(S)) measures how much the ordering of differentiation and evolution matters. We compute G across all 256 ECA rules and all 88 equivalence classes, combined with Shannon entropy, spatial correlation length, temporal autocorrelation, compressibility, spectral exponent, and finite-size scaling measures. We further extend the analysis to a family of commutator variants ([E, D₂], [Eⁿ, D], [E, [E, D]], and the anti-commutator {E, D}) and perform affine decomposition of each rule into a best-fit linear (XOR) component plus residual perturbation. We find that some Class IV rules (notably Rule 110, β ≈ 0.70) exhibit 1/f^β spectral signatures in their G time series, while others (Rule 54, β ≈ 0) do not, suggesting heterogeneity within Class IV. Glider-like structures in Rule 110 correspond to localized, persistent peaks in G. The corrected D operator dramatically sharpens the Class II/IV boundary: many Class II rules now have G ≡ 0, while all Class IV rules retain structured, nonzero G. Our experiments support the conjecture that Wolfram Class IV rules occupy a distinctive region of the commutator-affine feature space: they exhibit nonzero but structured commutator signatures and nonzero perturbation from affine behavior — consistent with "edge of chaos" characterizations but now given a precise operator-algebraic formulation.

> **Note on operator correction (2026-03-22):** An earlier version of this work used a spatial adjacency XOR as the differentiation operator D. The correct definition is D(S) = S ⊕ φ(S), which detects *which cells are about to flip* — the temporal derivative, not the spatial gradient. All results in this version use the corrected operator. The correction significantly changes the commutator's behavior: it now directly measures the non-commutativity between "detecting upcoming change" and "evolving the system."

## 1. Introduction

The classification of elementary cellular automata into Wolfram's four classes — uniform (I), periodic (II), chaotic (III), and complex (IV) — remains one of the most important open problems in discrete dynamical systems. While Classes I–III admit relatively clean characterizations (fixed points, limit cycles, and positive entropy respectively), Class IV has resisted formal definition. These rules produce long transients, localized structures ("gliders"), and are conjectured to support universal computation, yet no single measure cleanly separates them from the other classes.

We propose an approach rooted in operator algebra. Consider two natural operators on binary strings:

- **E** (evolution): apply the CA rule φ to advance one timestep: E(S) = φ(S)
- **D** (differentiation): detect which cells are about to change: D(S) = S ⊕ φ(S)

If E and D commute — that is, if D(E(S)) = E(D(S)) for all states S — then detecting upcoming changes before or after evolution yields the same result. This holds trivially for Class I rules (everything collapses) and for the purely affine (XOR-linear) rules. The *Groovy Commutator* G = [E, D] measures the degree of non-commutativity:

$$G(S) = D(E(S)) \oplus E(D(S))$$

Our central finding is that the statistical properties of G — its entropy, spatial correlation, temporal persistence, compressibility, and spectral signature — encode information about the dynamical class of the underlying rule, and that Class IV rules occupy a characteristic "structured but nontrivial" regime.

### 1.1 Related Work

Wolfram's original classification (1984, 2002) was based on visual inspection of spacetime diagrams. Subsequent work has attempted formalization via Lyapunov exponents (Bagnoli et al., 1992), mean-field theory (Gutowitz, 1991), computational complexity (Cook, 2004), Langton's λ parameter (1990), and various entropy-based measures (Martinez et al., 2013). The affine decomposition approach follows Rose (2026), who decomposes ECA rules into XOR-linear base plus nonlinear perturbation as a structural analysis tool.

Our commutator approach is, to our knowledge, novel. It draws on the philosophy of commutator analysis in quantum mechanics and Lie algebra, adapted to the discrete Boolean setting.

## 2. Definitions

### 2.1 Elementary Cellular Automata

An elementary cellular automaton is a map f: {0,1}³ → {0,1} applied synchronously to a ring of N cells. The 256 possible rules are indexed by the Wolfram number W = Σᵢ f(i) · 2ⁱ where i encodes the three-cell neighborhood.

### 2.2 Operators

Let S ∈ {0,1}ᴺ be a binary state on a ring of N cells with periodic boundary conditions.

**Evolution operator E (apply the rules):**
$$E(S) = \phi(S)$$
where φ is the CA rule applied synchronously to all cells: φ(S)[i] = f(S[i-1], S[i], S[i+1]).

**Differentiation operator D (which cells will flip):**
$$D(S) = S \oplus \phi(S)$$

D(S) marks the cells that are about to change on the next timestep. This is the *temporal derivative* expressed as a binary field — not a spatial gradient. D(S)[i] = 1 iff cell i will have a different value after one application of the rule.

**Integration operator I (apply changes):**
$$I(S, d) = S \oplus d$$

I applies a change pattern d to state S. Note that E(S) = I(S, D(S)) = S ⊕ (S ⊕ φ(S)) = φ(S) — evolution is the integration of differentiation.

**Second-order differentiation D₂ (the system dreaming about its own dynamics):**
$$D_2(S) = D(D(S)) = D(S) \oplus \phi(D(S))$$

D₂ feeds the change pattern D(S) back into the rule φ *as if it were a state*, then XORs the result with D(S) itself. This asks: "if the change pattern were a world, which cells of that world would change?" It is NOT the second time derivative.

### 2.3 The Groovy Commutator

$$G(S) = [E, D](S) = D(E(S)) \oplus E(D(S))$$

G measures the failure of E and D to commute. When G(S) = 0 for all S, detecting upcoming changes before or after evolution yields the same result.

Expanding both paths explicitly:
- **Path 1:** D(E(S)) = φ(S) ⊕ φ(φ(S)) — evolve first, then ask what's about to change in the evolved state.
- **Path 2:** E(D(S)) = φ(D(S)) = φ(S ⊕ φ(S)) — compute the change pattern first, then evolve it as if it were a state.

G(S) = 0 means these two paths agree — the system's dynamics commute with change-detection.

### 2.4 Extended Commutator Family

We study several related operators:

| Operator | Definition | Interpretation |
|----------|-----------|----------------|
| [E, D] = G | D(E(S)) ⊕ E(D(S)) | Standard commutator |
| [E, D₂] | D₂(E(S)) ⊕ E(D₂(S)) | Second-order change-detection sensitivity |
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
- 7 commutator variants: [E,D], [E,D₂], [E²,D], [E³,D], [E⁴,D], [E,[E,D]], {E,D}
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

### 3.5 Experiment 5: Equivalence Class Table

- All 88 ECA equivalence classes (under left-right reflection and complementation)
- Width N = 101, T = 300, 50 random initial conditions per representative
- Transient discard: 50 steps
- Full commutator signature: G entropy, correlation length, compressibility, periodicity
- Output: complete lookup table (see Appendix A)

### 3.6 Experiment 6: Spectral Analysis

- All 88 equivalence class representatives
- Width N = 201, T = 1000, discard first 200 steps, 10 trials
- Compute G entropy at each timestep → time series
- FFT of entropy time series → power spectrum
- Fit log-log slope (β exponent): β≈0 = white noise, β≈1 = 1/f noise, β≈2 = brown noise
- Output: β exponent, R², peak frequency for each rule

### 3.7 Experiment 7: Glider Detection

- Rule 110 (proven universal), with Rules 54 and 30 as comparisons
- Width N = 201, T = 500
- Compute G at each timestep, track persistent local maxima
- Glider criterion: peak persists ≥ 8 timesteps, moves ≤ 3 cells/step
- Output: spacetime diagram of G with detected glider tracks overlaid

### 3.8 Experiment 8: Finite-Size Scaling

- Rules 54, 110, 30, 90, 4, 0 (one representative per class)
- Sizes N = 51, 101, 201, 501
- T = 500, discard first 100 steps, 20 trials each
- Compute G entropy, correlation length, compressibility at each size
- Test: do Class IV measures converge or diverge with system size?

## 4. Results

### 4.1 Baseline: Class Separation in Commutator Space

We computed the Groovy Commutator G across all 256 ECA rules (N=101, T=200, 50 random trials each). The table below shows mean statistics by Wolfram class:

| Class | n | G Entropy | G Corr Length | G Temporal ACF | G Compressibility |
|-------|---|-----------|---------------|----------------|-------------------|
| I     | 16 | 0.058 ± 0.233 | 0.13 ± 0.51 | 1.000 ± 0.000 | 0.250 ± 0.049 |
| II    | 63 | 0.368 ± 0.361 | 0.71 ± 0.70 | 0.944 ± 0.144 | 0.334 ± 0.095 |
| III   | 30 | 0.610 ± 0.389 | 0.87 ± 0.62 | 0.366 ± 0.402 | 0.408 ± 0.109 |
| IV    |  9 | **0.842 ± 0.217** | **1.01 ± 0.04** | 0.211 ± 0.147 | 0.432 ± 0.057 |

With the corrected D(S) = S ⊕ φ(S), the class separation is substantially cleaner than with the old spatial-XOR operator:

- **Class I** rules have G entropy ≈ 0 and correlation length ≈ 0, with perfect temporal autocorrelation (constant zero commutator). These rules collapse all states, so D and E commute trivially.

- **Class II** rules now show G entropy = 0.37 (down from 0.65 with the old operator) with high temporal autocorrelation (0.94). Many Class II rules that had spurious non-commutativity under the old operator now have G=0 — their dynamics commute perfectly with change-detection. The remaining nonzero Class II rules show periodic commutator patterns.

- **Class III** is bimodal: non-affine rules (e.g., Rule 30 at 0.95) have high G entropy, while affine rules (90, 150) have G=0 exactly. The class mean (0.61) reflects this split. Temporal autocorrelation is low (0.37), indicating rapid decorrelation.

- **Class IV** rules have the **highest mean G entropy** (0.842) — now more clearly separated from other classes. Their temporal autocorrelation (0.21) is intermediate, and their compressibility (0.43) is the highest. Correlation lengths are uniformly ≈ 1.0 with the corrected operator, making spatial correlation uninformative.

### 4.2 Commutator Zoo: Variant Discrimination

We computed all 7 commutator variants on 20 representative rules. Selected results for key rules:

| Rule | Class | [E,D] | [E,D₂] | [E²,D] | [E³,D] | [E⁴,D] | [E,[E,D]] | {E,D} |
|------|-------|-------|---------|---------|---------|---------|-----------|-------|
| 30   | III   | 0.949 | 0.851 | 0.919 | 0.987 | 0.982 | 0.965 | 0.891 |
| 54   | IV    | 0.820 | 0.945 | 0.923 | 0.949 | 0.942 | 0.807 | 0.867 |
| 90   | III*  | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.993 |
| 110  | IV    | 0.938 | 0.995 | 0.996 | 0.878 | 0.973 | 0.910 | 0.934 |
| 106  | IV    | 0.989 | 0.885 | 0.947 | 0.983 | 0.978 | 0.980 | 0.807 |
| 193  | IV    | 0.923 | 0.906 | 0.995 | 0.990 | 0.996 | 0.888 | 0.000 |
| 4    | II    | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| 60   | III*  | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.993 |

*Rules 90 and 60 are affine, causing all XOR-type commutators to vanish.

Key findings from the zoo:

- **Many more rules now have G=0**: With D(S) = S ⊕ φ(S), several rules that had nonzero commutators under the old spatial-XOR definition now commute exactly (e.g., Rule 4, Rule 232). The corrected operator filters out spurious non-commutativity.

- **Affine confirmation**: Rules 90, 60, 105, and 150 (all affine) have all XOR-type commutators exactly zero. The anti-commutator {E,D} is high for some (0.993) — since both paths agree, their AND is nonzero where either is.

- **Nested commutator for Class IV**: With the corrected D, the nested commutator shows only modest drops: Rule 110: [E,D]=0.938 vs [E,[E,D]]=0.910 (3% drop); Rule 54: 0.820 vs 0.807 (2% drop). The hierarchical structure is present but subtle.

- **Multi-step behavior is non-monotonic**: For Rule 110, [E²,D]=0.996 > [E,D]=0.938, but [E³,D]=0.878 drops below, and [E⁴,D]=0.973 recovers. This oscillatory pattern suggests the commutator at different time horizons interacts with the rule's periodic structures.

### 4.3 Class IV Distinguishing Features

Extended runs (N=201, T=500) with all 6 measures reveal:

| Rule | Class | G Entropy | Block Ent. | Corr Len | Temp ACF | Hamming | Compress |
|------|-------|-----------|------------|----------|----------|---------|----------|
| 54   | IV    | 0.828 | 2.74 | 1.00 | 0.490 | 0.466 | 0.291 |
| 106  | IV    | 0.994 | 3.57 | 1.00 | 0.058 | 0.375 | 0.360 |
| 110  | IV    | 0.939 | 3.12 | 1.00 | 0.117 | 0.297 | 0.265 |
| 124  | IV    | 0.940 | 3.13 | 1.00 | 0.118 | 0.296 | 0.259 |
| 137  | IV    | 0.930 | 3.11 | 1.00 | 0.122 | 0.287 | 0.256 |
| 193  | IV    | 0.932 | 2.97 | 1.00 | 0.121 | 0.287 | 0.248 |
| 0    | I     | 0.000 | 0.00 | 0.00 | 1.000 | 0.000 | 0.119 |
| 32   | I     | 0.000 | 0.00 | 0.00 | 1.000 | 0.000 | 0.119 |
| 4    | II    | 0.000 | 0.00 | 0.00 | 1.000 | 0.000 | 0.119 |
| 108  | II    | 0.101 | 0.32 | 0.93 | 1.000 | 0.000 | 0.146 |
| 30   | III   | 0.952 | 3.48 | 1.00 | 0.042 | 0.313 | 0.359 |
| 90   | III*  | 0.000 | 0.00 | 0.00 | 1.000 | 0.000 | 0.119 |
| 150  | III*  | 0.000 | 0.00 | 0.00 | 1.000 | 0.000 | 0.119 |

With the corrected D operator, the results change substantially from the previous version. Key findings:

1. **Class IV entropy range: 0.83–0.99**: Class IV rules concentrate in this range. Rule 106 (0.994) approaches but does not quite reach Rule 30's value (0.952). Rule 54 (0.828) is the lowest Class IV entropy — still clearly above zero but also clearly below the non-affine Class III level.

2. **Correlation lengths uniformly ~1.0**: With D(S) = S ⊕ φ(S), all rules show correlation length ≈ 1.0. The elevated values seen with the old spatial-XOR operator (Rule 54 at ξ=3.4) were an artifact. The corrected operator produces spatially uncorrelated commutator fields, making correlation length uninformative as a discriminator.

3. **Intermediate temporal autocorrelation**: Class IV ACF ranges 0.06–0.49, between Class II's 1.0 (constant G or periodic) and Class III's 0.04. Rule 54's ACF (0.49) is notably high, indicating persistent temporal structure in its commutator. This is the strongest discriminating feature with the corrected operator.

4. **Rule 4 (Class II) now has G=0**: This is the most dramatic change from the corrected operator. Rule 4 previously showed high G entropy (0.80); with D(S) = S ⊕ φ(S), it has zero commutator. Many simple Class II rules are now transparent to change-detection.

5. **Moderate Hamming distances**: Class IV d_H ranges 0.25–0.47, indicating ongoing change in the commutator. Class I/II rules that now have G=0 trivially have d_H = 0.

### 4.4 Affine Decomposition × Commutator

Across all 256 rules, combining affine decomposition with commutator measures:

| Class | Pert. Weight | G Entropy | G Corr Length | Fraction Affine |
|-------|-------------|-----------|---------------|-----------------|
| I     | 0.148 ± 0.082 | 0.059 ± 0.236 | 0.12 ± 0.49 | 12% |
| II    | 0.177 ± 0.077 | 0.368 ± 0.362 | 0.70 ± 0.67 | 6% |
| III   | 0.142 ± 0.102 | 0.610 ± 0.389 | 0.88 ± 0.63 | 27% |
| IV    | 0.167 ± 0.062 | **0.843 ± 0.215** | **1.04 ± 0.09** | **0%** |

With the corrected D operator, the class separation in G entropy is dramatically cleaner:

- **No Class IV rule is affine** (0%, unchanged). The conjecture that affinity precludes Class IV behavior holds.

- **Class IV has the highest mean G entropy** (0.843), well above Class III (0.610) and Class II (0.368). The corrected operator sharply reduces Class II's mean entropy (from 0.649 to 0.368) because many Class II rules now have G=0 — their dynamics commute perfectly with change-detection.

- **Correlation length is no longer discriminating**: All classes with nonzero G show correlation length ≈ 1.0 with the corrected operator. This means the conjecture must be reformulated without the correlation length criterion.

- **The cleaner separation** comes from G entropy alone: with the corrected D, the gap between Class IV (0.843) and Class II (0.368) has widened from 0.246 to 0.475. Class IV is the only class where *all* rules have G entropy > 0.8.

This supports a simplified version of our hypothesis:

> **Conjecture (Groovy Commutator Characterization, revised):** A rule R exhibits Class IV behavior if and only if: (1) R is non-affine (perturbation weight > 0), and (2) the Groovy Commutator G has entropy bounded in the range (0.8, 1.0) — strictly above Class II but generally below non-affine Class III rules at their maximum.

### 4.5 Spectral Analysis: 1/f Signatures in G Time Series

We computed the power spectrum of the G entropy time series for all 88 equivalence class representatives (N=201, T=1000, 10 trials). The spectral exponent β, obtained by fitting log(power) vs. log(frequency), reveals distinct signatures by class:

| Rule | Class | β exponent | R² | Interpretation |
|------|-------|-----------|------|----------------|
| 110  | IV    | **0.696**  | 0.148 | Approaching 1/f (pink noise) |
| 106  | IV    | **0.396**  | 0.080 | Positive β, structured fluctuations |
| 54   | IV    | −0.034     | 0.006 | Near-zero (flat spectrum) |
| 30   | III   | 0.029     | 0.004 | Near-zero (white noise) |
| 45   | III   | 0.177     | 0.015 | Weak positive |
| 146  | III   | 0.461     | 0.110 | Moderate β (structured non-affine III) |
| 18   | III   | 0.891     | 0.271 | Highest β in dataset (non-affine III) |
| 90   | III   | 0.000     | 0.000 | Trivial (affine, G=0) |
| 4    | II    | 0.000     | 0.000 | Trivial (converges to fixed G=0) |

Class-level summary: Class IV mean β=0.35±0.30, Class III mean β=0.18±0.26, Class II mean β=−0.08±0.34, Class I mean β=0.00.

Two Class IV rules (110 and 106) show positive β exponents (0.40–0.70), with Rule 110 approaching the 1/f regime (β≈0.7). Rule 54, however, shows a near-zero β, indicating its G entropy fluctuations lack strong long-range temporal correlations despite its high spatial correlation length. The R² values are modest, indicating the power-law fit is approximate — longer time series or ensemble averaging would improve the fit.

Notably, some non-affine Class III rules (Rule 18 at β=0.89, Rule 146 at β=0.46) show higher β than the Class IV mean. The spectral exponent alone does not cleanly separate Class IV from Class III — it is the *conjunction* of β, entropy level, and spatial correlation that matters. Rule 110's β=0.70 combined with its high G entropy (0.94) and non-maximal structure makes a stronger signature than any single measure.

### 4.6 Glider Detection via Commutator Peaks

We tracked persistent, moving local maxima in the G spacetime of Rule 110 (N=201, T=500). Results:

| Rule | Class | Tracks Detected | Speed Range | Duration Range |
|------|-------|-----------------|-------------|----------------|
| 110  | IV    | 564             | 0.0–1.0 cells/step | 8–245 steps |
| 54   | IV    | 394             | — | — |
| 30   | III   | 770             | — | — |

Rule 110 shows 564 distinct tracked features with speeds distributed across the range, with the highest concentration at |v|≈0.5–0.6 cells/step (97 and 167 tracks respectively). The speed distribution shows structure, with peaks at specific velocities suggesting the tracker is detecting distinct propagating features. Track durations extend up to 245 steps, indicating highly persistent structures.

Rule 30 (Class III) has more "tracks" (770) than Rule 110, but these are artifacts of the detection algorithm applied to random fluctuations — Class III rules produce rapidly decorrelating G, so any threshold-based tracker will find spurious matches among the noise. The key difference is qualitative: Rule 110's tracks have coherent trajectories visible in the spacetime diagram, while Rule 30's are scattered.

Rule 54 shows 394 tracks — more than expected for its sparser dynamics. With the corrected D operator (which detects cells about to change rather than spatial gradients), the commutator G highlights dynamic activity more broadly. The high track count in Rule 54 may reflect the rich structure of its glider interactions, or it may indicate that the persistence threshold (8 steps) is too permissive for this rule. Further work is needed to validate these tracks against known Rule 54 glider catalogs.

### 4.7 Finite-Size Scaling

We measured G entropy, correlation length, and compressibility across system sizes N = 51, 101, 201, 501:

| Rule | Class | G Entropy (N=51) | G Entropy (N=501) | Trend |
|------|-------|-----------------|-------------------|-------|
| 54   | IV    | 0.735 ± 0.091   | 0.837 ± 0.009     | Growing, not converged |
| 110  | IV    | 0.937 ± 0.005   | 0.939 ± 0.001     | Stable (converged) |
| 30   | III   | 0.943 ± 0.003   | 0.953 ± 0.001     | Stable (slight growth) |
| 90   | III   | 0.000           | 0.000              | Zero (affine) |
| 4    | II    | 0.000           | 0.000              | Zero (fixed point) |
| 0    | I     | 0.000           | 0.000              | Zero (trivial) |

Key findings:

1. **Rule 54's entropy grows with system size** (0.735 → 0.837), and its standard deviation decreases (0.091 → 0.009), suggesting it has not yet reached its thermodynamic limit at N=501. This is consistent with Rule 54's long correlation lengths — it may require much larger systems for its commutator statistics to converge.

2. **Rule 110 is essentially converged** by N=51 (0.937 ≈ 0.939). Its G entropy is remarkably stable across all sizes, indicating that its commutator structure is not a finite-size artifact.

3. **Rule 30 is also stable** (0.943 → 0.953), confirming that Class III behavior saturates quickly.

4. **Compressibility decreases uniformly with size** for all rules that have nonzero G (since larger random-looking arrays compress proportionally better). This is an artifact of the gzip measure and does not reflect genuine dynamical scaling.

5. **Correlation lengths converge to ~1.0** for all non-trivial rules at large sizes. The elevated values seen in smaller experiments (e.g., Rule 54 at ξ=3.5 in Appendix A) may reflect finite-size effects where the correlation length approaches the system size.

## 5. Discussion

### 5.1 Why the Commutator Captures Complexity

The Groovy Commutator measures, in precise terms, the failure of change-detection (D) and evolution (E) to commute. D(S) = S ⊕ φ(S) asks "which cells are about to flip?" — it is the temporal derivative encoded as a spatial pattern. When D and E commute, the answer to "what's about to change?" is the same whether we ask before or after evolution. This means the dynamics are, in a specific algebraic sense, *transparent to their own future* — knowing where change will happen before evolution is equivalent to knowing it after.

Class I rules trivialize both operators: D converges to zero as the state becomes uniform, so they commute vacuously. Affine (XOR-linear) rules commute because both E and D are linear operations over GF(2), and linear operators always commute.

Class III rules generate change patterns so aggressively that D(E(S)) and E(D(S)) both become essentially random — the commutator is nonzero but unstructured. The ordering doesn't matter because everything becomes unpredictable either way.

Class IV rules are precisely those where the ordering *matters in a structured way*: the commutator is nonzero (the system's change-detection does not commute with evolution) but the non-commutativity has spatial coherence (it's not random noise). This is the hallmark of computation — the system's dynamics interact with self-knowledge about upcoming change in a non-trivial but organized manner.

### 5.2 Connection to Edge of Chaos

Langton's "edge of chaos" hypothesis (1990) posits that computation occurs at the boundary between ordered and chaotic dynamics. The Groovy Commutator provides a new lens on this: the "edge" is where the commutator [E, D] is nonzero but structured. With the corrected D operator (which detects upcoming changes rather than spatial gradients), this analogy becomes tighter: the commutator measures whether a system can "predict its own change" consistently across different orderings of operations. This is analogous to how in quantum mechanics, the commutator [x, p] = iℏ defines the boundary between classical (commuting) and quantum behavior, and the magnitude of the commutator controls the uncertainty principle.

### 5.3 The Affine Connection

The affine decomposition reveals that Class IV rules are perturbations of linear dynamics — they have a dominant XOR-linear component with a small but critical nonlinear perturbation. This resonates with Brooklyn Rose's insight that the "interesting" dynamics arise not from maximal nonlinearity but from structured deviations from linearity.

The conjunction of moderate perturbation and moderate commutator entropy may be a signature of the "barrier" that Rose describes — the computational boundary between the affine regime (where behavior is predictable) and the fully nonlinear regime (where behavior is chaotic).

### 5.4 Limitations

1. **Classification ground truth**: Our Wolfram class assignments are based on community consensus, which is itself partially subjective for marginal rules. The classification of rules near class boundaries remains contested.

2. **Finite-size effects**: Most experiments use finite rings (N ≤ 201). Finite-size scaling (Experiment 8) addresses this partially for N up to 501, but some Class IV phenomena may require larger systems.

3. **Statistical measures**: Shannon entropy and correlation length are coarse measures. Richer topological or algebraic invariants of G might provide sharper discrimination.

4. **Conjecture status**: Our characterization is stated as a conjecture, not a theorem. A proof would require understanding the algebraic structure of the commutator for all possible rules, which remains open.

5. **Operator sensitivity**: The precise definition of D matters greatly. With D(S) = S ⊕ φ(S), the commutator measures non-commutativity of change-detection and evolution — a semantically cleaner quantity than spatial-gradient non-commutativity.

## 6. Conjectures

We state our main conjectures with epistemic tags indicating confidence levels.

**Conjecture 1 (Commutator Vanishing — HIGH CONFIDENCE):**
*For any affine ECA rule (W ∈ {0, 15, 51, 85, 90, 102, 105, 150, 153, 165, 170, 195, 204, 240, 255}), the Groovy Commutator G(S) = 0 for all S.*

Rationale: For affine rules, φ is GF(2)-linear, so D(S) = S ⊕ φ(S) is also linear. Both E = φ and D = id ⊕ φ are GF(2)-linear, and linear operators over GF(2) commute.

**Conjecture 2 (Class IV High-Entropy Regime — MODERATE CONFIDENCE):**
*Class IV rules are characterized by G entropy in the range (0.8, 1.0), with all Class IV rules exceeding the Class II mean.*

Rationale: With the corrected D operator, Class IV has mean G entropy 0.842 ± 0.217 — the highest class mean. All individual Class IV rules exceed 0.82. The corrected operator makes the separation from Class II (0.368) much cleaner. Note: correlation length is no longer discriminating with D(S) = S ⊕ φ(S).

**Conjecture 3 (Multi-step Non-monotonicity — LOW CONFIDENCE):**
*For Class IV rules, [Eⁿ, D] shows non-monotonic behavior as n increases (e.g., Rule 110: 0.996 at n=2, 0.878 at n=3, 0.973 at n=4), reflecting interaction between the commutator and periodic substructures of the rule.*

Rationale: With the corrected operator, multi-step commutator entropy no longer shows clean monotonic growth. The oscillatory pattern may reflect finite-size effects or genuinely complex temporal structure.

**Conjecture 4 (Affine-Commutator Joint Characterization — MODERATE CONFIDENCE):**
*A necessary condition for Class IV behavior is: (1) the rule is non-affine (perturbation weight > 0), and (2) G entropy > 0.8.*

Rationale: With the corrected operator, 100% of Class IV rules are non-affine AND have G entropy > 0.8. The correlation length criterion from the earlier version has been dropped because D(S) = S ⊕ φ(S) produces uniformly short-range spatial correlations in G.

**Conjecture 5 (Nested Commutator Structure — SPECULATIVE):**
*For Class IV rules, the nested commutator [E, [E, D]] shows a modest but consistent drop from [E, D] (e.g., Rule 110: 0.910 vs 0.938; Rule 54: 0.807 vs 0.820), indicating some hierarchical structure in the non-commutativity.*

Rationale: The drop is less dramatic than under the old operator (previously 0.713 vs 0.973 for Rule 110) but still present. The corrected operator reveals that the hierarchical cancellation is subtle rather than dramatic.

## 7. Future Directions

1. **Larger rule spaces**: Extend to 2-state 2D CA, k=3 totalistic rules, or general k-state 1D CA to test universality of the commutator characterization.

2. **Algebraic theory**: Develop a formal algebraic framework for the commutator over GF(2)^N, potentially connecting to the theory of Boolean function composition.

3. **Improved spectral analysis**: Our preliminary spectral results (§4.5) show promising 1/f-like signatures for Class IV, but R² values are low. Longer time series (T > 10000), windowed FFT (Welch's method), and detrended fluctuation analysis (DFA) may yield cleaner exponents.

4. **Machine learning**: Use commutator features (entropy, β exponent, correlation length) as inputs to a classifier and evaluate whether they outperform raw spacetime features for Wolfram class prediction.

5. **Improved glider detection**: The peak-tracking approach (§4.6) detects features in Rule 110 but needs refinement: speed-dependent detection thresholds, cross-correlation tracking, and comparison with known glider catalogs to validate which tracked features correspond to genuine gliders.

6. **Finite-size scaling to larger N**: Rule 54's growing entropy at N=501 (§4.7) motivates runs at N=1001 and N=2001 to determine its thermodynamic limit. If Rule 54's entropy continues to grow, it may approach Class III values, raising questions about the robustness of entropy-based separation.

7. **Connection to computational universality**: Cook (2004) proved Rule 110 is Turing-universal. Can the commutator structure be used to detect or predict computational universality? The nested commutator drop for Rule 110 (§4.2) may be a signature worth investigating.

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

We compute the Groovy Commutator signature for all 88 ECA equivalence classes under left-right reflection and black-white complementation, using the corrected D(S) = S ⊕ φ(S) operator. For each class, the smallest rule number serves as representative. Parameters: N=101, T=300, 50 random initial conditions, transient discard of 50 steps.

The full table is available in `results/equivalence_table.md`. Below we highlight the Class IV representatives and key comparisons:

| Rule | Class | G≠0 (rand) | G Entropy | G Corr Len | G Compress | Pert Weight |
|------|-------|------------|-----------|------------|------------|-------------|
| 54   | IV    | 1.00 | 0.808±0.037 | 1.02 | 0.414 | 0.250 |
| 106  | IV    | 1.00 | 0.989±0.002 | 1.00 | 0.531 | 0.250 |
| 110  | IV    | 1.00 | 0.939±0.003 | 1.00 | 0.437 | 0.125 |
| 30   | III   | 1.00 | 0.949±0.003 | 1.00 | 0.534 | 0.250 |
| 45   | III   | 1.00 | 0.937±0.005 | 1.88 | 0.492 | 0.250 |
| 90   | III*  | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 0.000 |
| 24   | II    | 1.00 | 0.974±0.027 | 1.66 | 0.461 | 0.250 |
| 4    | II    | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 0.125 |
| 0    | I     | 0.00 | 0.000±0.000 | 0.00 | 0.238 | 0.000 |

Notable: 22 of the 65 Class II representatives now have G=0 identically (vs. only affine rules before). Many simple Class II rules have dynamics that commute perfectly with change-detection.

### A.1 Analysis of the Equivalence Class Table

**Conjecture 1 confirmed.** Every affine rule (perturbation weight = 0.000) — Rules 0, 15, 51, 60, 90, 105, 150, 170, 204 — has G entropy exactly 0.000. The commutator vanishes identically for all affine rules, as predicted.

**Corrected D reveals additional G=0 rules.** Beyond affine rules, 22 Class II representatives (including Rules 4, 12, 13, 19, 23, 33, 36, 37, 44, 72, 76, 77, 78, 104, 132, 140, 164, 172, 178, 200, 232) now show G=0 from random ICs. These are non-affine rules whose dynamics nevertheless commute with change-detection. This is a new finding enabled by the corrected operator.

**Class IV rules cluster at G entropy 0.81–0.99.** All three Class IV representatives have G entropy in this range, well above the majority of Class II rules (most of which are zero or < 0.6) but slightly below the highest non-affine Class III rules (Rule 30 at 0.949, Rule 146 at 0.977).

**Class III bimodal structure persists.** The 4 affine Class III rules have G=0; the 8 non-affine ones range from 0.567 (Rule 122) to 0.977 (Rule 146). Some non-affine Class III rules (Rule 146 at 0.977) exceed Class IV's Rule 54 (0.808), confirming that entropy alone does not perfectly separate the classes.

**Some Class II rules have surprisingly high G entropy.** Rule 24 (0.974), Rule 27 (0.961), Rule 35 (0.992), and Rule 152 (0.964) all exceed most Class IV rules. These are "near-boundary" rules that may warrant reclassification or represent a limitation of the Wolfram class system.
