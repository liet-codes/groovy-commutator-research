"""
Experiment 6: Spectral Analysis of G Time Series

For each of the 88 equivalence class representatives:
- Run N=201, T=1000, discard first 200 steps
- Compute G at each timestep, record Shannon entropy of G over time
- Compute power spectrum (FFT) of the G entropy time series
- Fit log-log slope (beta exponent) to look for 1/f^beta signatures
  - beta ~ 0 = white noise (Class III?)
  - beta ~ 1 = pink/1/f noise (Class IV?)
  - beta ~ 2 = brown noise (Class II?)

Output: results/spectral.csv, results/spectral_examples.png
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from src.eca import step, random_initial
from src.operators import G as groovy_commutator
from src.measures import shannon_entropy
from experiments.exp5_equivalence_table import get_equivalence_classes, WOLFRAM_CLASS


def compute_spectral_exponent(time_series):
    """Compute the spectral exponent beta from a time series.

    Fits log(power) vs log(frequency) to get the slope beta.
    Returns (beta, r_squared, peak_frequency).
    """
    n = len(time_series)
    if n < 16:
        return 0.0, 0.0, 0.0

    # Remove mean
    ts = time_series - np.mean(time_series)

    # Check if signal is essentially constant
    if np.std(ts) < 1e-10:
        return 0.0, 0.0, 0.0

    # Compute power spectrum via FFT
    fft_vals = np.fft.rfft(ts)
    power = np.abs(fft_vals) ** 2
    freqs = np.fft.rfftfreq(n)

    # Skip DC component (index 0) and use positive frequencies
    freqs = freqs[1:]
    power = power[1:]

    # Remove zero-power entries for log
    mask = power > 0
    if np.sum(mask) < 5:
        return 0.0, 0.0, 0.0

    log_freq = np.log10(freqs[mask])
    log_power = np.log10(power[mask])

    # Fit linear regression in log-log space
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_freq, log_power)
    beta = -slope  # Convention: power ~ 1/f^beta, so slope is negative

    # Peak frequency (excluding DC)
    peak_idx = np.argmax(power)
    peak_freq = freqs[peak_idx]

    return beta, r_value ** 2, peak_freq


def run():
    N = 201
    T = 1000
    T_DISCARD = 200
    N_TRIALS = 10
    SEED = 42

    classes = get_equivalence_classes()
    reps = sorted(classes.keys())
    print(f"Experiment 6: Spectral Analysis of G Time Series")
    print(f"N={N}, T={T}, discard={T_DISCARD}, trials={N_TRIALS}")
    print(f"Analyzing {len(reps)} equivalence class representatives")

    rng = np.random.default_rng(SEED)
    results = []

    for rep in reps:
        wclass = WOLFRAM_CLASS.get(rep, '?')

        trial_betas = []
        trial_r2s = []
        trial_peaks = []

        for trial in range(N_TRIALS):
            state = random_initial(N, rng)

            # Evolve and collect G entropy time series
            g_entropy_series = []
            for t in range(T):
                g = groovy_commutator(state, rep)
                if t >= T_DISCARD:
                    g_entropy_series.append(shannon_entropy(g))
                state = step(state, rep)

            g_entropy_series = np.array(g_entropy_series)

            beta, r2, peak_freq = compute_spectral_exponent(g_entropy_series)
            trial_betas.append(beta)
            trial_r2s.append(r2)
            trial_peaks.append(peak_freq)

        mean_beta = np.mean(trial_betas)
        mean_r2 = np.mean(trial_r2s)
        mean_peak = np.mean(trial_peaks)

        results.append({
            'rule': rep,
            'class': wclass,
            'beta_exponent': round(mean_beta, 4),
            'r_squared': round(mean_r2, 4),
            'peak_frequency': round(mean_peak, 6),
        })
        print(f"  Rule {rep:>3} (Class {wclass:>3}): "
              f"beta={mean_beta:.3f} R²={mean_r2:.3f} peak_f={mean_peak:.4f}")

    # Write CSV
    os.makedirs('results', exist_ok=True)
    with open('results/spectral.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['rule', 'class', 'beta_exponent', 'r_squared', 'peak_frequency'])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved results/spectral.csv")

    # Generate example plots: one rule per class
    example_rules = {
        'I': 0,      # Class I
        'II': 4,     # Class II
        'III': 30,   # Class III (non-affine)
        'IV': 110,   # Class IV
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    class_order = ['I', 'II', 'III', 'IV']

    for ax, cls_name in zip(axes.flat, class_order):
        rule = example_rules[cls_name]
        state = random_initial(N, np.random.default_rng(42))

        g_entropy_series = []
        for t in range(T):
            g = groovy_commutator(state, rule)
            if t >= T_DISCARD:
                g_entropy_series.append(shannon_entropy(g))
            state = step(state, rule)

        ts = np.array(g_entropy_series)
        ts_centered = ts - np.mean(ts)

        if np.std(ts_centered) < 1e-10:
            ax.text(0.5, 0.5, f'Rule {rule} (Class {cls_name})\nConstant signal (G=0)',
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
            ax.set_title(f'Rule {rule} — Class {cls_name}')
            continue

        fft_vals = np.fft.rfft(ts_centered)
        power = np.abs(fft_vals) ** 2
        freqs = np.fft.rfftfreq(len(ts))

        # Skip DC
        freqs = freqs[1:]
        power = power[1:]

        mask = power > 0
        log_freq = np.log10(freqs[mask])
        log_power = np.log10(power[mask])

        ax.scatter(log_freq, log_power, s=2, alpha=0.5, color='steelblue')

        # Fit line
        slope, intercept, r_val, _, _ = stats.linregress(log_freq, log_power)
        fit_line = slope * log_freq + intercept
        beta = -slope
        ax.plot(log_freq, fit_line, 'r-', linewidth=2,
                label=f'β = {beta:.2f} (R² = {r_val**2:.3f})')

        ax.set_xlabel('log₁₀(frequency)')
        ax.set_ylabel('log₁₀(power)')
        ax.set_title(f'Rule {rule} — Class {cls_name}')
        ax.legend(fontsize=10)

    plt.suptitle('Power Spectra of G Entropy Time Series', fontsize=14)
    plt.tight_layout()
    plt.savefig('results/spectral_examples.png', dpi=150)
    print("Saved results/spectral_examples.png")
    plt.close()

    # Summary by class
    print("\n=== Summary by Wolfram Class ===")
    for cls in ['I', 'II', 'III', 'IV']:
        cls_results = [r for r in results if r['class'] == cls]
        if cls_results:
            betas = [r['beta_exponent'] for r in cls_results]
            r2s = [r['r_squared'] for r in cls_results]
            print(f"  Class {cls:>3} (n={len(cls_results):>2}): "
                  f"beta={np.mean(betas):.3f}±{np.std(betas):.3f}  "
                  f"R²={np.mean(r2s):.3f}")


if __name__ == '__main__':
    run()
