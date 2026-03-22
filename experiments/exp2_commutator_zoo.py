"""Experiment 2: Compare all commutator variants on interesting rules."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.eca import evolve, random_initial
from src.operators import COMMUTATOR_VARIANTS
from src.measures import shannon_entropy, spatial_correlation_length

# Interesting rules spanning all classes
INTERESTING_RULES = [
    0, 4, 22, 30, 32, 45, 54, 60, 73, 90, 105, 106, 108, 110,
    122, 126, 150, 184, 193, 232
]

WIDTH = 101
TIMESTEPS = 200
N_TRIALS = 30
TRANSIENT = 50


def run_zoo():
    print("Running Experiment 2: Commutator Zoo")
    print(f"Rules: {INTERESTING_RULES}")
    print(f"Variants: {list(COMMUTATOR_VARIANTS.keys())}")

    rng = np.random.default_rng(42)
    records = []

    for rule in INTERESTING_RULES:
        for variant_name, variant_fn in COMMUTATOR_VARIANTS.items():
            entropies = []
            corr_lens = []

            for trial in range(N_TRIALS):
                init = random_initial(WIDTH, rng)
                history = evolve(init, rule, TIMESTEPS)

                trial_entropies = []
                trial_corrs = []
                for t in range(TRANSIENT, TIMESTEPS):
                    c_state = variant_fn(history[t], rule)
                    trial_entropies.append(shannon_entropy(c_state))
                    trial_corrs.append(spatial_correlation_length(c_state))

                entropies.append(np.mean(trial_entropies))
                corr_lens.append(np.mean(trial_corrs))

            records.append({
                'rule': rule,
                'variant': variant_name,
                'entropy_mean': np.mean(entropies),
                'entropy_std': np.std(entropies),
                'corr_length_mean': np.mean(corr_lens),
                'corr_length_std': np.std(corr_lens),
            })

        print(f"  Completed rule {rule}")

    df = pd.DataFrame(records)
    os.makedirs('results', exist_ok=True)
    df.to_csv('results/zoo.csv', index=False)
    print(f"Saved results/zoo.csv ({len(df)} records)")

    # Generate heatmap: rules × variants, colored by entropy
    pivot = df.pivot(index='rule', columns='variant', values='entropy_mean')
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    im = ax.imshow(pivot.values, aspect='auto', cmap='viridis')
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=9)
    ax.set_xlabel('Commutator Variant')
    ax.set_ylabel('ECA Rule')
    ax.set_title('Commutator Zoo: Entropy by Rule × Variant')
    plt.colorbar(im, ax=ax, label='Entropy')
    plt.tight_layout()
    plt.savefig('results/zoo_heatmap.png', dpi=150)
    print("Saved results/zoo_heatmap.png")
    plt.close()

    return df


if __name__ == '__main__':
    run_zoo()
