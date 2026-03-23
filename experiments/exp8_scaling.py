"""
Experiment 8: Finite-Size Scaling of Commutator Measures

Run key rules at multiple system sizes to check whether Class IV
measures converge or diverge with system size.

Rules: 54, 110, 30, 90, 4, 0
Sizes: N = 51, 101, 201, 501
T = 500, 20 trials each

Output: results/scaling.csv, results/scaling_plot.png
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.eca import step, random_initial
from src.operators import G as groovy_commutator
from src.measures import shannon_entropy, spatial_correlation_length, compressibility_ratio

RULES = [54, 110, 30, 90, 4, 0]
SIZES = [51, 101, 201, 501]
T = 500
T_DISCARD = 100
N_TRIALS = 20
SEED = 42

RULE_CLASS = {
    54: 'IV', 110: 'IV',
    30: 'III', 90: 'III',
    4: 'II', 0: 'I',
}


def run():
    print(f"Experiment 8: Finite-Size Scaling")
    print(f"Rules: {RULES}")
    print(f"Sizes: {SIZES}")
    print(f"T={T}, discard={T_DISCARD}, trials={N_TRIALS}")

    rng = np.random.default_rng(SEED)
    results = []

    for rule in RULES:
        cls = RULE_CLASS[rule]
        for N in SIZES:
            trial_entropy = []
            trial_corr = []
            trial_compress = []

            for trial in range(N_TRIALS):
                state = random_initial(N, rng)

                g_entropies = []
                g_corrs = []
                g_compresses = []

                for t in range(T):
                    g = groovy_commutator(state, rule)
                    if t >= T_DISCARD:
                        g_entropies.append(shannon_entropy(g))
                        g_corrs.append(spatial_correlation_length(g))
                        g_compresses.append(compressibility_ratio(g))
                    state = step(state, rule)

                trial_entropy.append(np.mean(g_entropies))
                trial_corr.append(np.mean(g_corrs))
                trial_compress.append(np.mean(g_compresses))

            result = {
                'rule': rule,
                'class': cls,
                'N': N,
                'g_entropy_mean': round(np.mean(trial_entropy), 6),
                'g_entropy_std': round(np.std(trial_entropy), 6),
                'g_corr_mean': round(np.mean(trial_corr), 4),
                'g_corr_std': round(np.std(trial_corr), 4),
                'g_compress_mean': round(np.mean(trial_compress), 6),
                'g_compress_std': round(np.std(trial_compress), 6),
            }
            results.append(result)
            print(f"  Rule {rule:>3} (Class {cls:>3}) N={N:>3}: "
                  f"ent={result['g_entropy_mean']:.4f}±{result['g_entropy_std']:.4f} "
                  f"corr={result['g_corr_mean']:.2f}±{result['g_corr_std']:.2f} "
                  f"comp={result['g_compress_mean']:.4f}")

    # Write CSV
    os.makedirs('results', exist_ok=True)
    with open('results/scaling.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved results/scaling.csv")

    # Generate scaling plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    measures = [
        ('g_entropy_mean', 'g_entropy_std', 'G Entropy'),
        ('g_corr_mean', 'g_corr_std', 'G Correlation Length'),
        ('g_compress_mean', 'g_compress_std', 'G Compressibility'),
    ]

    class_colors = {'I': 'blue', 'II': 'green', 'III': 'red', 'IV': 'gold'}
    markers = {0: 's', 4: '^', 30: 'D', 54: 'o', 90: 'v', 110: '*'}

    for ax, (mean_col, std_col, title) in zip(axes, measures):
        for rule in RULES:
            cls = RULE_CLASS[rule]
            rule_data = [r for r in results if r['rule'] == rule]
            sizes = [r['N'] for r in rule_data]
            means = [r[mean_col] for r in rule_data]
            stds = [r[std_col] for r in rule_data]

            ax.errorbar(sizes, means, yerr=stds,
                        marker=markers[rule], color=class_colors[cls],
                        label=f'R{rule} ({cls})', linewidth=2, markersize=8,
                        capsize=4)

        ax.set_xlabel('System Size N')
        ax.set_ylabel(title)
        ax.set_title(title)
        ax.set_xscale('log')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Finite-Size Scaling of Commutator Measures', fontsize=14)
    plt.tight_layout()
    plt.savefig('results/scaling_plot.png', dpi=150)
    print("Saved results/scaling_plot.png")
    plt.close()

    # Summary
    print("\n=== Scaling Summary ===")
    for rule in RULES:
        cls = RULE_CLASS[rule]
        rule_data = [r for r in results if r['rule'] == rule]
        ent_small = [r['g_entropy_mean'] for r in rule_data if r['N'] == 51][0]
        ent_large = [r['g_entropy_mean'] for r in rule_data if r['N'] == 501][0]
        corr_small = [r['g_corr_mean'] for r in rule_data if r['N'] == 51][0]
        corr_large = [r['g_corr_mean'] for r in rule_data if r['N'] == 501][0]
        print(f"  Rule {rule:>3} ({cls:>3}): "
              f"entropy {ent_small:.4f}→{ent_large:.4f} "
              f"corr_len {corr_small:.2f}→{corr_large:.2f}")


if __name__ == '__main__':
    run()
