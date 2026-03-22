"""
Experiment 5: Equivalence Class Lookup Table
Compute commutator signatures for all 88 ECA equivalence classes.

For each class:
- G ≠ 0 from single black cell
- G ≠ 0 fraction from random ICs
- G entropy (mean ± std)
- G correlation length
- G compressibility
- G periodicity (detected period or "aperiodic")
- Affine perturbation magnitude
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import csv
from src.eca import step
from src.operators import G as groovy_commutator_op
from src.measures import shannon_entropy, spatial_correlation_length, compressibility_ratio as compressibility
from src.affine import affine_decompose

# 88 equivalence class representatives (standard ordering)
# Under left-right reflection R and complement C, each rule has up to 4 equivalents
# We pick the smallest rule number as representative

def get_equivalence_classes():
    """Return dict mapping representative -> list of equivalent rules."""
    seen = set()
    classes = {}
    for rule in range(256):
        if rule in seen:
            continue
        # Compute equivalents
        equiv = set()
        # Original
        equiv.add(rule)
        # Left-right reflection
        reflected = reflect_rule(rule)
        equiv.add(reflected)
        # Complement
        complemented = complement_rule(rule)
        equiv.add(complemented)
        # Both
        equiv.add(complement_rule(reflected))
        
        rep = min(equiv)
        classes[rep] = sorted(equiv)
        seen.update(equiv)
    return classes

def reflect_rule(rule):
    """Reflect a rule (swap left and right neighbors)."""
    new_rule = 0
    for i in range(8):
        L = (i >> 2) & 1
        C_val = (i >> 1) & 1
        R = i & 1
        reflected_i = (R << 2) | (C_val << 1) | L
        bit = (rule >> i) & 1
        new_rule |= (bit << reflected_i)
    return new_rule

def complement_rule(rule):
    """Complement a rule (swap 0 and 1)."""
    new_rule = 0
    for i in range(8):
        complement_i = 7 - i
        bit = (rule >> i) & 1
        new_rule |= ((1 - bit) << complement_i)
    return new_rule

# Known Wolfram classifications for representatives
WOLFRAM_CLASS = {
    0: 'I', 1: 'II', 2: 'II', 3: 'II', 4: 'II', 5: 'II',
    6: 'II', 7: 'II', 8: 'I', 9: 'II', 10: 'II', 11: 'II',
    12: 'II', 13: 'II', 14: 'II', 15: 'II', 18: 'III', 19: 'II',
    22: 'III', 23: 'II', 24: 'II', 25: 'II', 26: 'II', 27: 'II',
    28: 'II', 29: 'II', 30: 'III', 32: 'I', 33: 'II', 34: 'II',
    35: 'II', 36: 'II', 37: 'II', 38: 'II', 40: 'I', 41: 'II',
    42: 'II', 43: 'II', 44: 'II', 45: 'III', 46: 'II', 50: 'II',
    51: 'II', 54: 'IV', 56: 'II', 57: 'II', 58: 'II', 60: 'III',
    62: 'II', 72: 'II', 73: 'II', 74: 'II', 76: 'II', 77: 'II',
    78: 'II', 90: 'III', 94: 'II', 104: 'II', 105: 'III',
    106: 'IV', 108: 'II', 110: 'IV', 122: 'III', 126: 'III',
    128: 'I', 130: 'III', 132: 'II', 134: 'II', 136: 'I',
    138: 'II', 140: 'II', 142: 'II', 146: 'III', 150: 'III',
    152: 'II', 154: 'II', 156: 'II', 160: 'I', 162: 'II',
    164: 'II', 168: 'I', 170: 'II', 172: 'II', 178: 'II',
    184: 'II', 200: 'II', 204: 'II', 232: 'II',
}

def detect_periodicity(series, max_period=50):
    """Detect periodicity in a binary series. Returns period or 0 for aperiodic."""
    if len(series) < 10:
        return 0
    # Check if all zero
    if all(s == 0 for s in series):
        return 1  # trivially periodic
    
    for p in range(1, min(max_period, len(series)//3)):
        is_periodic = True
        for i in range(p, min(3*p, len(series))):
            if series[i] != series[i % p]:
                is_periodic = False
                break
        if is_periodic:
            return p
    return 0  # aperiodic

def run():
    N = 101
    T = 300
    T_DISCARD = 50
    N_RANDOM = 50

    classes = get_equivalence_classes()
    print(f"Found {len(classes)} equivalence classes")

    results = []

    for rep in sorted(classes.keys()):
        equiv = classes[rep]
        wclass = WOLFRAM_CLASS.get(rep, '?')

        # --- Single black cell test ---
        single_cell = np.zeros(N, dtype=np.uint8)
        single_cell[N // 2] = 1
        state = single_cell.copy()

        g_nonzero_single = False
        for t in range(T):
            g = groovy_commutator_op(state, rep)
            if np.any(g):
                g_nonzero_single = True
                break
            state = step(state, rep)

        # --- Random IC tests ---
        g_nonzero_count = 0
        entropies = []
        corr_lengths = []
        compressibilities = []
        hamming_series_entropies = []
        periodicity_results = []

        for trial in range(N_RANDOM):
            state = np.random.randint(0, 2, N, dtype=np.uint8)
            
            g_ever_nonzero = False
            g_entropy_series = []
            
            for t in range(T):
                g = groovy_commutator_op(state, rep)
                
                if t >= T_DISCARD:
                    if np.any(g):
                        g_ever_nonzero = True
                    g_entropy_series.append(shannon_entropy(g))
                
                state = step(state, rep)
            
            if g_ever_nonzero:
                g_nonzero_count += 1
            
            if g_entropy_series:
                entropies.append(np.mean(g_entropy_series))
                
                # Compute final G for correlation/compression
                final_state = state
                final_g = groovy_commutator_op(final_state, rep)
                corr_lengths.append(spatial_correlation_length(final_g))
                compressibilities.append(compressibility(final_g))
                
                # Periodicity: quantize entropy series
                quantized = [1 if e > 0.5 else 0 for e in g_entropy_series]
                period = detect_periodicity(quantized)
                periodicity_results.append(period)

        # Affine decomposition
        ad = affine_decompose(rep)
        affine_base = ad['best_affine_name']
        pert_weight = ad['perturbation_weight']

        result = {
            'representative': rep,
            'equivalents': '|'.join(map(str, equiv)),
            'wolfram_class': wclass,
            'n_equivalents': len(equiv),
            'g_nonzero_single_cell': g_nonzero_single,
            'g_nonzero_random_frac': g_nonzero_count / N_RANDOM if N_RANDOM > 0 else 0,
            'g_entropy_mean': np.mean(entropies) if entropies else 0,
            'g_entropy_std': np.std(entropies) if entropies else 0,
            'g_corr_length_mean': np.mean(corr_lengths) if corr_lengths else 0,
            'g_compressibility_mean': np.mean(compressibilities) if compressibilities else 0,
            'g_periodicity_mode': max(set(periodicity_results), key=periodicity_results.count) if periodicity_results else 0,
            'affine_base': affine_base,
            'perturbation_weight': pert_weight,
        }
        results.append(result)
        
        period_str = str(result['g_periodicity_mode']) if result['g_periodicity_mode'] > 0 else 'aperiodic'
        print(f"  Rule {rep:>3} (Class {wclass:>3}): "
              f"G_ent={result['g_entropy_mean']:.3f} "
              f"G_corr={result['g_corr_length_mean']:.2f} "
              f"G_single={'Y' if g_nonzero_single else 'N'} "
              f"G_rand={result['g_nonzero_random_frac']:.2f} "
              f"pert={pert_weight:.3f} "
              f"period={period_str}")

    # Write CSV
    os.makedirs('results', exist_ok=True)
    with open('results/equivalence_classes.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    # Write markdown table for paper
    with open('results/equivalence_table.md', 'w') as f:
        f.write("# ECA Equivalence Class Commutator Signatures\n\n")
        f.write("| Rep | Class | Equiv | G≠0 (single) | G≠0 (random) | G Entropy | G Corr Len | G Compress | Period | Pert Weight |\n")
        f.write("|-----|-------|-------|---------------|--------------|-----------|------------|------------|--------|-------------|\n")
        for r in results:
            period_str = str(r['g_periodicity_mode']) if r['g_periodicity_mode'] > 0 else 'aper.'
            f.write(f"| {r['representative']:>3} | {r['wolfram_class']:>3} "
                    f"| {r['n_equivalents']} "
                    f"| {'✓' if r['g_nonzero_single_cell'] else '✗'} "
                    f"| {r['g_nonzero_random_frac']:.2f} "
                    f"| {r['g_entropy_mean']:.3f}±{r['g_entropy_std']:.3f} "
                    f"| {r['g_corr_length_mean']:.2f} "
                    f"| {r['g_compressibility_mean']:.3f} "
                    f"| {period_str} "
                    f"| {r['perturbation_weight']:.3f} |\n")

    print(f"\nResults written to results/equivalence_classes.csv")
    print(f"Markdown table written to results/equivalence_table.md")

if __name__ == '__main__':
    run()
