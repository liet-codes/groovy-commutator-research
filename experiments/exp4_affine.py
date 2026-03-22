"""Experiment 4: Affine decomposition × commutator analysis."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.eca import evolve, random_initial
from src.operators import G
from src.measures import shannon_entropy, spatial_correlation_length
from src.affine import affine_decompose, perturbation_coherence

# Wolfram class assignments (same as exp1)
WOLFRAM_CLASS = {}
CLASS_I = [0, 8, 32, 40, 128, 136, 160, 168, 255, 223, 239, 247, 254, 250, 252, 253]
CLASS_II = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 19, 23, 24, 25, 26, 27, 28, 29,
            33, 34, 35, 36, 37, 38, 42, 43, 44, 46, 50, 51, 56, 57, 58, 62, 72, 76, 77, 78,
            94, 104, 108, 130, 132, 134, 138, 140, 142, 152, 154, 156, 162, 164, 170, 172,
            178, 184, 200, 204, 232]
CLASS_III = [18, 22, 30, 45, 60, 73, 75, 86, 89, 90, 101, 102, 105, 107, 109, 121, 122,
             126, 129, 146, 149, 150, 153, 161, 165, 169, 181, 182, 183, 195]
CLASS_IV = [41, 54, 106, 110, 124, 131, 137, 147, 193]

for r in CLASS_I: WOLFRAM_CLASS[r] = 1
for r in CLASS_II: WOLFRAM_CLASS[r] = 2
for r in CLASS_III: WOLFRAM_CLASS[r] = 3
for r in CLASS_IV: WOLFRAM_CLASS[r] = 4

WIDTH = 101
TIMESTEPS = 200
N_TRIALS = 30
TRANSIENT = 50


def run_affine():
    print("Running Experiment 4: Affine Decomposition × Commutator")
    rng = np.random.default_rng(42)
    records = []

    for rule in range(256):
        # Affine decomposition
        dec = affine_decompose(rule)
        p_coherence = perturbation_coherence(dec['perturbation'])

        # Commutator measures
        entropies = []
        corr_lens = []

        for trial in range(N_TRIALS):
            init = random_initial(WIDTH, rng)
            history = evolve(init, rule, TIMESTEPS)

            g_ent = []
            for t in range(TRANSIENT, TIMESTEPS):
                g_state = G(history[t], rule)
                g_ent.append(shannon_entropy(g_state))

            entropies.append(np.mean(g_ent))
            corr_lens.append(spatial_correlation_length(G(history[TIMESTEPS], rule)))

        records.append({
            'rule': rule,
            'wolfram_class': WOLFRAM_CLASS.get(rule, 0),
            'best_affine': dec['best_affine_name'],
            'perturbation_magnitude': dec['perturbation_magnitude'],
            'perturbation_weight': dec['perturbation_weight'],
            'perturbation_coherence': p_coherence,
            'is_affine': dec['is_affine'],
            'g_entropy': np.mean(entropies),
            'g_entropy_std': np.std(entropies),
            'g_corr_length': np.mean(corr_lens),
        })

        if rule % 32 == 0:
            print(f"  Completed rule {rule}/255")

    df = pd.DataFrame(records)
    os.makedirs('results', exist_ok=True)
    df.to_csv('results/affine.csv', index=False)
    print(f"Saved results/affine.csv ({len(df)} rules)")

    # Generate 2D scatter: perturbation magnitude vs G entropy, colored by class
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    class_colors = {0: 'gray', 1: 'blue', 2: 'green', 3: 'red', 4: 'gold'}
    class_labels = {0: 'Unclassified', 1: 'Class I', 2: 'Class II', 3: 'Class III', 4: 'Class IV'}

    # Plot 1: Perturbation magnitude vs G entropy
    ax = axes[0]
    for cls in [0, 1, 2, 3, 4]:
        mask = df['wolfram_class'] == cls
        if mask.any():
            ax.scatter(df.loc[mask, 'perturbation_weight'],
                       df.loc[mask, 'g_entropy'],
                       c=class_colors[cls], label=class_labels[cls],
                       alpha=0.7, s=30, edgecolors='k', linewidths=0.3)
    ax.set_xlabel('Perturbation Weight (distance from affine)')
    ax.set_ylabel('G Entropy')
    ax.set_title('Affine Perturbation vs Groovy Commutator')
    ax.legend(fontsize=8)

    # Plot 2: Perturbation coherence vs G entropy
    ax = axes[1]
    for cls in [0, 1, 2, 3, 4]:
        mask = df['wolfram_class'] == cls
        if mask.any():
            ax.scatter(df.loc[mask, 'perturbation_coherence'],
                       df.loc[mask, 'g_entropy'],
                       c=class_colors[cls], label=class_labels[cls],
                       alpha=0.7, s=30, edgecolors='k', linewidths=0.3)
    ax.set_xlabel('Perturbation Coherence')
    ax.set_ylabel('G Entropy')
    ax.set_title('Perturbation Coherence vs Groovy Commutator')
    ax.legend(fontsize=8)

    plt.suptitle('Hypothesis Test: Class IV = Nonzero Perturbation AND Structured Commutator', fontsize=12)
    plt.tight_layout()
    plt.savefig('results/affine_scatter.png', dpi=150)
    print("Saved results/affine_scatter.png")
    plt.close()

    # Print hypothesis test results
    print("\n=== HYPOTHESIS TEST ===")
    print("Class IV = nonzero perturbation AND nonzero structured commutator?")
    for cls in [1, 2, 3, 4]:
        mask = df['wolfram_class'] == cls
        if mask.any():
            sub = df[mask]
            print(f"\nClass {cls}:")
            print(f"  Perturbation weight: {sub['perturbation_weight'].mean():.3f} ± {sub['perturbation_weight'].std():.3f}")
            print(f"  G entropy: {sub['g_entropy'].mean():.4f} ± {sub['g_entropy'].std():.4f}")
            print(f"  G corr length: {sub['g_corr_length'].mean():.2f} ± {sub['g_corr_length'].std():.2f}")
            print(f"  Fraction affine: {sub['is_affine'].mean():.2f}")

    return df


if __name__ == '__main__':
    run_affine()
