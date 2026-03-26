"""Experiment 10: Block CA Commutator vs Coarse-Grained Effective Rule

Tests whether the Groovy Commutator detects structure at the sub-coarse-graining
scale that the effective-rule collapse erases.

Key test: Block Rule 30 collapses to ECA Rule 150 (affine, G≡0).
Does Block Rule 30 on Margolus topology show nonzero G?

Reference: Brooklyn Rose, "Elementary Interaction Graphs" (2026), §5-6.
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.operators import G as eca_G, phi
from src.measures import shannon_entropy


# ═══ Vectorized Block CA Engine ═══

def block_tables(block_rule: int):
    """Return L and R truth tables for a block rule."""
    R_num = block_rule % 16
    L_num = block_rule // 16
    L = np.array([(L_num >> i) & 1 for i in range(4)], dtype=np.uint8)
    R = np.array([(R_num >> i) & 1 for i in range(4)], dtype=np.uint8)
    return L, R


def block_step(state, L, R, phase):
    """Vectorized block step. phase=0 for even, phase=1 for odd."""
    n = len(state)
    result = state.copy()
    start = phase
    # Vectorized: compute all block indices at once
    positions = np.arange(start, n - 1, 2)
    a = state[positions]
    b = state[positions + 1]
    idx = a * 2 + b
    result[positions] = L[idx]
    result[positions + 1] = R[idx]
    return result


def block_phi(state, L, R):
    """One full Margolus step: even + odd phase."""
    s = block_step(state, L, R, 0)
    return block_step(s, L, R, 1)


def block_G(state, L, R):
    """G(S) for block CA."""
    e_s = block_phi(state, L, R)
    e_e_s = block_phi(e_s, L, R)
    d_of_e = e_s ^ e_e_s
    d_s = state ^ e_s
    e_of_d = block_phi(d_s, L, R)
    return d_of_e ^ e_of_d


def measure_block_g(block_rule, n_samples=50, N=100, T=200, skip=50, seed=42):
    """Measure G entropy for a block rule."""
    rng = np.random.default_rng(seed)
    L, R = block_tables(block_rule)
    entropies = []
    densities = []
    
    for _ in range(n_samples):
        state = rng.integers(0, 2, size=N, dtype=np.uint8)
        g_ents = []
        g_dens = []
        for t in range(T):
            g = block_G(state, L, R)
            g_ents.append(shannon_entropy(g))
            g_dens.append(np.mean(g))
            state = block_phi(state, L, R)
        entropies.append(np.mean(g_ents[skip:]))
        densities.append(np.mean(g_dens[skip:]))
    
    return np.mean(entropies), np.std(entropies), np.mean(densities)


def measure_eca_g(eca_rule, n_samples=50, N=100, T=200, skip=50, seed=42):
    """Measure G entropy for an ECA rule."""
    rng = np.random.default_rng(seed)
    entropies = []
    
    for _ in range(n_samples):
        state = rng.integers(0, 2, size=N, dtype=np.uint8)
        g_ents = []
        for t in range(T):
            g = eca_G(state, eca_rule)
            g_ents.append(shannon_entropy(g))
            state = phi(state, eca_rule)
        entropies.append(np.mean(g_ents[skip:]))
    
    return np.mean(entropies), np.std(entropies)


def main():
    print("=" * 70)
    print("EXPERIMENT 10: Block CA Commutator vs Effective Rule")
    print("=" * 70)
    
    # Verify Block Rule 30 = rule 198
    L, R = block_tables(198)
    assert list(L) == [0, 0, 1, 1], f"L wrong: {L}"  # L(x,y) = x
    assert list(R) == [0, 1, 1, 0], f"R wrong: {R}"  # R(x,y) = x⊕y
    print("Block Rule 30 (rule 198): L(x,y)=x, R(x,y)=x⊕y ✓")
    print("Effective rule: a⊕b⊕c = ECA Rule 150 (affine, Class III, G≡0)")
    
    # ═══ Key Test: Block Rule 30 vs ECA Rule 150 ═══
    print("\n--- KEY TEST: Block Rule 30 vs ECA Rule 150 ---")
    
    print("  Running Block Rule 30 (Margolus topology)...")
    b_ent, b_std, b_dens = measure_block_g(198)
    print(f"  Block Rule 30 G entropy: {b_ent:.6f} ± {b_std:.6f} (density: {b_dens:.4f})")
    
    print("  Running ECA Rule 150 (elementary topology)...")
    e_ent, e_std = measure_eca_g(150)
    print(f"  ECA Rule 150 G entropy:  {e_ent:.6f} ± {e_std:.6f}")
    
    if b_ent > 0.01 and e_ent < 0.01:
        print(f"\n  ★ BLOCK RULE 30 HAS NONZERO G ({b_ent:.4f}), ECA RULE 150 HAS G≡0")
        print("    → GC detects phase-alternation structure invisible to effective rule!")
    elif b_ent < 0.01:
        print(f"\n  Both have G≈0. GC consistent with effective-rule collapse.")
    else:
        print(f"\n  Both nonzero. Difference: {abs(b_ent - e_ent):.6f}")
    
    # ═══ Control tests ═══
    print("\n--- Control: Additional block rules ---")
    
    controls = [
        (202, "Block Identity (L=x, R=y)"),
        (172, "Block Swap (L=y, R=x)"),
        (136, "Block AND (L=x∧y, R=x∧y)"),
        (0,   "Block Zero (L=0, R=0)"),
        (255, "Block All-Ones (L=1, R=1)"),
    ]
    
    for br, label in controls:
        ent, std, dens = measure_block_g(br, n_samples=20, T=100)
        print(f"  {label} (rule {br}): G entropy = {ent:.6f} ± {std:.6f}")
    
    # ═══ Scan all 256 block rules (faster with vectorized engine) ═══
    print("\n--- Scanning all 256 block rules ---")
    print("  (20 samples × 100 steps each)")
    
    results = []
    for br in range(256):
        ent, std, dens = measure_block_g(br, n_samples=10, N=80, T=100, skip=20)
        results.append((br, ent))
    
    nonzero = [(r, g) for r, g in results if g > 0.01]
    zero = [(r, g) for r, g in results if g <= 0.01]
    
    print(f"  Block rules with G≈0: {len(zero)}/256")
    print(f"  Block rules with G>0: {len(nonzero)}/256")
    
    if nonzero:
        sorted_nz = sorted(nonzero, key=lambda x: -x[1])
        print(f"\n  Top 15 block rules by G entropy:")
        for r, g in sorted_nz[:15]:
            L_num, R_num = r // 16, r % 16
            print(f"    Block rule {r:>3} (L={L_num:>2}, R={R_num:>2}): G = {g:.4f}")
    
    # ═══ Summary ═══
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Block Rule 30 (Margolus) G entropy: {b_ent:.6f}")
    print(f"  ECA Rule 150 (elementary) G entropy: {e_ent:.6f}")
    print(f"  Block rules with G>0: {len(nonzero)}/256 ({100*len(nonzero)/256:.1f}%)")
    print(f"  Block rules with G≈0: {len(zero)}/256 ({100*len(zero)/256:.1f}%)")
    
    if b_ent > 0.01:
        print(f"\n  ★ The Groovy Commutator detects structure at the Margolus phase-")
        print(f"    alternation level that is invisible to the effective-rule collapse.")
        print(f"    The GC is FINER than compositional-level equivalence.")


if __name__ == '__main__':
    main()
