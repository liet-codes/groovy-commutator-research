"""Experiment 3: Focused Class IV analysis with longer runs."""

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
    temporal_autocorrelation_mean, hamming_distance_normalized,
    compressibility_ratio, block_entropy
)

# Rules to compare
CLASS_IV_RULES = [54, 106, 110, 124, 137, 193]
CLASS_I_RULES = [0, 32]
CLASS_II_RULES = [4, 108]
CLASS_III_RULES = [30, 90, 150]
ALL_RULES = CLASS_IV_RULES + CLASS_I_RULES + CLASS_II_RULES + CLASS_III_RULES

RULE_CLASS = {}
for r in CLASS_IV_RULES: RULE_CLASS[r] = 4
for r in CLASS_I_RULES: RULE_CLASS[r] = 1
for r in CLASS_II_RULES: RULE_CLASS[r] = 2
for r in CLASS_III_RULES: RULE_CLASS[r] = 3

WIDTH = 201
TIMESTEPS = 500
N_TRIALS = 30
TRANSIENT = 100


def run_class4_hunt():
    print("Running Experiment 3: Class IV Hunt")
    print(f"Width={WIDTH}, Timesteps={TIMESTEPS}, Trials={N_TRIALS}")
    print(f"Rules: {ALL_RULES}")

    rng = np.random.default_rng(42)
    records = []

    for rule in ALL_RULES:
        cls = RULE_CLASS[rule]
        all_g_entropy = []
        all_g_block_entropy = []
        all_corr_len = []
        all_temporal_acf = []
        all_hamming = []
        all_compress = []

        for trial in range(N_TRIALS):
            init = random_initial(WIDTH, rng)
            history = evolve(init, rule, TIMESTEPS)

            g_ent_series = []
            g_states = []
            for t in range(TRANSIENT, TIMESTEPS):
                g_state = G(history[t], rule)
                g_states.append(g_state)
                g_ent_series.append(shannon_entropy(g_state))

            all_g_entropy.append(np.mean(g_ent_series))

            # Block entropy of final G state
            if g_states:
                all_g_block_entropy.append(block_entropy(g_states[-1], block_size=4))

            # Correlation length
            if g_states:
                all_corr_len.append(spatial_correlation_length(g_states[-1]))

            # Temporal autocorrelation
            if len(g_ent_series) > 10:
                all_temporal_acf.append(temporal_autocorrelation_mean(np.array(g_ent_series)))

            # Hamming distance between successive G states
            if len(g_states) > 1:
                hd = [hamming_distance_normalized(g_states[i], g_states[i + 1])
                      for i in range(len(g_states) - 1)]
                all_hamming.append(np.mean(hd))

            # Compressibility
            if g_states:
                all_compress.append(compressibility_ratio(g_states[-1]))

        records.append({
            'rule': rule,
            'wolfram_class': cls,
            'g_entropy': np.mean(all_g_entropy),
            'g_entropy_std': np.std(all_g_entropy),
            'g_block_entropy': np.mean(all_g_block_entropy) if all_g_block_entropy else 0,
            'g_corr_length': np.mean(all_corr_len) if all_corr_len else 0,
            'g_temporal_acf': np.mean(all_temporal_acf) if all_temporal_acf else 0,
            'g_hamming': np.mean(all_hamming) if all_hamming else 0,
            'g_compressibility': np.mean(all_compress) if all_compress else 0,
        })
        print(f"  Rule {rule} (Class {cls}): entropy={records[-1]['g_entropy']:.4f}")

    df = pd.DataFrame(records)
    os.makedirs('results', exist_ok=True)
    df.to_csv('results/class4_hunt.csv', index=False)
    print(f"Saved results/class4_hunt.csv")

    # Generate comparison plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    measures = [
        ('g_entropy', 'G Entropy'),
        ('g_block_entropy', 'G Block Entropy (k=4)'),
        ('g_corr_length', 'G Correlation Length'),
        ('g_temporal_acf', 'G Temporal Autocorrelation'),
        ('g_hamming', 'G Hamming Distance (successive)'),
        ('g_compressibility', 'G Compressibility'),
    ]
    class_colors = {1: 'blue', 2: 'green', 3: 'red', 4: 'gold'}
    class_names = {1: 'I', 2: 'II', 3: 'III', 4: 'IV'}

    for ax, (col, title) in zip(axes.flat, measures):
        for cls in [1, 2, 3, 4]:
            mask = df['wolfram_class'] == cls
            vals = df.loc[mask, col].values
            rules = df.loc[mask, 'rule'].values
            ax.bar([f"R{r}" for r in rules], vals,
                   color=class_colors[cls], alpha=0.7,
                   label=f'Class {class_names[cls]}')
        ax.set_title(title, fontsize=10)
        ax.tick_params(axis='x', rotation=45, labelsize=7)
        ax.legend(fontsize=7)

    plt.suptitle('Class IV Hunt: Comparing Wolfram Classes', fontsize=14)
    plt.tight_layout()
    plt.savefig('results/class4_hunt_comparison.png', dpi=150)
    print("Saved results/class4_hunt_comparison.png")
    plt.close()

    return df


if __name__ == '__main__':
    run_class4_hunt()
