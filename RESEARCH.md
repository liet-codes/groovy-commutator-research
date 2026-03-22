# Groovy Commutator Research

## Goal
Formal characterization of the Groovy Commutator and its relationship to Wolfram's Class IV cellular automata. We're hunting for a **computational definition of Class IV** using commutator structure.

## Core Questions
1. Does the Groovy Commutator reliably distinguish Class IV from Classes I/II/III?
2. Are there other commutator-like functions (from the Boolean function space) that better characterize Class IV?
3. What is the spectral signature of Class IV rules under commutator analysis?
4. Can we define Class IV as a specific regime of commutator behavior?

## Definitions

### The Groovy Commutator
For a CA with evolution operator E and differentiation operator D:
```
G(S) = C( D(E(S)), E(D(S)) )
```
where C is a comparison function (typically XOR or difference).

In ECA context:
- **E**: Apply the CA rule to advance one timestep
- **D**: Compute the spatial derivative (XOR of adjacent cells, or difference)
- **C**: Compare two states (XOR, Hamming distance, etc.)

### Wolfram Classes (informal)
- **Class I**: Evolve to uniform state (G → 0)
- **Class II**: Evolve to periodic/oscillating patterns (G periodic)
- **Class III**: Chaotic, pseudo-random (G noisy, high entropy)
- **Class IV**: Complex, long transients, localized structures (G structured but non-periodic)

## Research Plan

### Phase 1: Baseline Characterization
- Implement all 256 ECA rules
- Compute G for each rule across many initial conditions
- Classify G behavior: zero, periodic, chaotic, structured
- Compare against known Wolfram classifications

### Phase 2: Commutator Zoo
- Systematically explore other commutator-like operators:
  - [E, D] = E∘D - D∘E (standard)
  - [E, D²] = second-order commutator
  - [E, I] where I is integration (cumulative XOR)
  - [E^n, D] for multi-step evolution
  - Nested: [E, [E, D]] (Jacobi-like)
  - Anti-commutator: {E, D} = E∘D + D∘E
- For each, measure:
  - Entropy of the commutator output
  - Spatial correlation length
  - Temporal autocorrelation
  - Lyapunov-like sensitivity (perturbation response)

### Phase 3: Spectral Analysis
- Compute power spectra of G time series
- Look for spectral signatures that distinguish Class IV
- Connect to graph Laplacian eigenvalues where possible

### Phase 4: Boolean Function Space
- ECA rules ARE Boolean functions of 3 inputs
- Decompose each rule into affine + non-linear perturbation (Brooklyn's approach)
- Map the commutator behavior against the affine decomposition
- Look for the conjunction: nonzero perturbation AND nonzero commutator AND coherence

### Phase 5: Paper
- Formal definitions
- Experimental results
- Conjectures
- Connection to existing literature

## Key References
- Wolfram, S. (2002). A New Kind of Science.
- Rose, B. (2026). The Shape of the Barrier.
- Rose, B. (2026). Affine Decompositions of Elementary Cellular Automata. (Brooklyn's paper)
- Cook, M. (2004). Universality in Elementary Cellular Automata. Complex Systems.
- Langton, C. (1990). Computation at the Edge of Chaos.
- Martinez, G. et al. (2013). On the Classification of Complex Behaviors in CAs.
