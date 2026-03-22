"""Experiment 1: Baseline characterization of all 256 ECA rules."""

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
from src.measures import (
    shannon_entropy, spatial_correlation_length,
    temporal_autocorrelation_mean, compressibility_ratio
)

# Known Wolfram classifications (representative assignments)
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
N_TRIALS = 50
TRANSIENT = 50  # skip initial transient


def run_baseline():
    print("Running Experiment 1: Baseline characterization of all 256 ECA rules")
    print(f"Width={WIDTH}, Timesteps={TIMESTEPS}, Trials={N_TRIALS}")

    rng = np.random.default_rng(42)
    records = []

    for rule in range(256):
        entropies = []
        corr_lens = []
        temporal_acfs = []
        compressibilities = []

        for trial in range(N_TRIALS):
            init = random_initial(WIDTH, rng)
            history = evolve(init, rule, TIMESTEPS)

            # Compute G at each timestep after transient
            g_entropies = []
            g_states = []
            for t in range(TRANSIENT, TIMESTEPS):
                g_state = G(history[t], rule)
                g_states.append(g_state)
                g_entropies.append(shannon_entropy(g_state))

            # Average entropy
            mean_entropy = np.mean(g_entropies)
            entropies.append(mean_entropy)

            # Correlation length of final G state
            if len(g_states) > 0:
                corr_lens.append(spatial_correlation_length(g_states[-1]))

            # Temporal autocorrelation of G entropy series
            if len(g_entropies) > 10:
                temporal_acfs.append(temporal_autocorrelation_mean(np.array(g_entropies)))

            # Compressibility of final G state
            if len(g_states) > 0:
                compressibilities.append(compressibility_ratio(g_states[-1]))

        wolfram_class = WOLFRAM_CLASS.get(rule, 0)
        records.append({
            'rule': rule,
            'wolfram_class': wolfram_class,
            'g_entropy_mean': np.mean(entropies),
            'g_entropy_std': np.std(entropies),
            'g_corr_length_mean': np.mean(corr_lens) if corr_lens else 0,
            'g_temporal_acf_mean': np.mean(temporal_acfs) if temporal_acfs else 0,
            'g_compressibility_mean': np.mean(compressibilities) if compressibilities else 0,
        })

        if rule % 32 == 0:
            print(f"  Completed rule {rule}/255")

    df = pd.DataFrame(records)
    os.makedirs('results', exist_ok=True)
    df.to_csv('results/baseline.csv', index=False)
    print(f"Saved results/baseline.csv ({len(df)} rules)")

    # Generate scatter plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    class_colors = {0: 'gray', 1: 'blue', 2: 'green', 3: 'red', 4: 'gold'}
    class_labels = {0: 'Unclassified', 1: 'Class I', 2: 'Class II', 3: 'Class III', 4: 'Class IV'}

    for cls in [0, 1, 2, 3, 4]:
        mask = df['wolfram_class'] == cls
        if mask.any():
            ax.scatter(df.loc[mask, 'g_entropy_mean'],
                       df.loc[mask, 'g_corr_length_mean'],
                       c=class_colors[cls], label=class_labels[cls],
                       alpha=0.7, s=30, edgecolors='k', linewidths=0.3)

    ax.set_xlabel('G Entropy (mean)')
    ax.set_ylabel('G Correlation Length (mean)')
    ax.set_title('Groovy Commutator: Entropy vs Correlation Length by Wolfram Class')
    ax.legend()
    plt.tight_layout()
    plt.savefig('results/baseline_scatter.png', dpi=150)
    print("Saved results/baseline_scatter.png")
    plt.close()

    return df


if __name__ == '__main__':
    run_baseline()
