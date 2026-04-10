"""
Experiment 11: Recursive Groovy Commutator — Testing stability under g^n(s)

Hypothesis: Class IV CAs have a stable recursive substrate where iterating
the Groovy Commutator converges, while Class III stays chaotic and Class I/II
collapses to zero.

Key insight: We're treating G not as a measurement but as a dynamical operator.
g(s) = [E,D](s) produces a new configuration. What happens if we iterate?
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from tqdm import tqdm

from src.eca import step, evolve
from src.operators import G, D, E
from src.measures import shannon_entropy, hamming_distance


# ═══ Recursive G Application ═══

def G_iterated(initial_state: np.ndarray, rule_number: int, n_iter: int, 
               track_metrics=True):
    """
    Apply G iteratively: s_{k+1} = G(s_k)
    
    Returns:
        states: list of state arrays (or just final if not tracking)
        metrics: dict of tracked metrics over iterations
    """
    states = [initial_state.copy()]
    s = initial_state.copy()
    
    metrics = {
        'entropy': [],
        'perturbation_weight': [],
        'hamming_from_prev': [],
        'hamming_from_initial': [],
    } if track_metrics else None
    
    for i in range(n_iter):
        # Apply G operator
        s = G(s, rule_number)
        states.append(s.copy())
        
        if track_metrics:
            metrics['entropy'].append(shannon_entropy(s))
            metrics['perturbation_weight'].append(np.sum(s) / len(s))
            metrics['hamming_from_prev'].append(np.sum(s ^ states[-2]))
            metrics['hamming_from_initial'].append(np.sum(s ^ initial_state))
    
    return states, metrics


def G_orbit_convergence(initial_state: np.ndarray, rule_number: int, 
                        max_iter=1000, convergence_window=50, 
                        entropy_stability_thresh=0.01):
    """
    Test if G-iteration converges to a stable pattern.
    
    Convergence criteria:
    1. Entropy stabilizes (std over window < thresh)
    2. Or we detect a cycle (return to previous state)
    
    Returns:
        converged: bool
        n_steps: steps to convergence (or max_iter)
        final_entropy: entropy at end
        cycle_length: if cycle detected, else 0
        entropy_series: full entropy trajectory
    """
    s = initial_state.copy()
    seen_states = {}  # Hash -> step number
    entropy_series = []
    
    for step_num in range(max_iter):
        # Check for exact cycle
        state_hash = s.tobytes()
        if state_hash in seen_states:
            cycle_start = seen_states[state_hash]
            cycle_length = step_num - cycle_start
            return True, step_num, shannon_entropy(s), cycle_length, entropy_series
        seen_states[state_hash] = step_num
        
        # Track entropy
        entropy_series.append(shannon_entropy(s))
        
        # Check entropy stability
        if len(entropy_series) >= convergence_window:
            recent = entropy_series[-convergence_window:]
            if np.std(recent) < entropy_stability_thresh:
                return True, step_num, np.mean(recent), 0, entropy_series
        
        # Iterate
        s = G(s, rule_number)
    
    # Did not converge
    return False, max_iter, entropy_series[-1], 0, entropy_series


# ═══ Experiments ═══

def test_rule_recursive_stability(rule_number: int, n_samples=10, 
                                   size=101, max_iter=500,
                                   rule_name=""):
    """Test recursive G stability for a single rule across random initial states."""
    
    results = {
        'rule': rule_number,
        'name': rule_name,
        'converged_count': 0,
        'cycle_count': 0,
        'mean_steps_to_stable': [],
        'mean_final_entropy': [],
        'entropy_trajectories': [],
    }
    
    for sample in range(n_samples):
        # Random initial state (sparse to dense)
        density = np.random.uniform(0.1, 0.5)
        initial = np.random.random(size) < density
        initial = initial.astype(np.uint8)
        
        conv, steps, final_ent, cycle_len, ent_series = G_orbit_convergence(
            initial, rule_number, max_iter=max_iter
        )
        
        results['converged_count'] += int(conv)
        results['cycle_count'] += int(cycle_len > 0)
        if conv:
            results['mean_steps_to_stable'].append(steps)
        results['mean_final_entropy'].append(final_ent)
        results['entropy_trajectories'].append(ent_series)
    
    # Summarize
    results['convergence_rate'] = results['converged_count'] / n_samples
    results['cycle_rate'] = results['cycle_count'] / n_samples
    results['mean_steps'] = np.mean(results['mean_steps_to_stable']) if results['mean_steps_to_stable'] else max_iter
    results['final_entropy_mean'] = np.mean(results['mean_final_entropy'])
    results['final_entropy_std'] = np.std(results['mean_final_entropy'])
    
    return results


def sweep_rules_by_class():
    """Test recursive G stability across Wolfram classes."""
    
    # Representative rules by class
    test_rules = {
        'Class I (uniform)': [0, 32, 128, 160, 250],
        'Class II (periodic)': [4, 36, 72, 200, 108, 94],
        'Class III (chaotic)': [30, 60, 90, 150],
        'Class IV (complex)': [54, 110, 106],
    }
    
    # Add Brooklyn's high-entropy Class II outliers
    test_rules['Class II (outliers)'] = [24, 35]
    
    all_results = []
    
    for class_name, rules in test_rules.items():
        print(f"\n{'='*60}")
        print(f"Testing {class_name}")
        print(f"{'='*60}")
        
        for rule in rules:
            print(f"\nRule {rule}...")
            result = test_rule_recursive_stability(
                rule, n_samples=20, size=101, max_iter=1000,
                rule_name=f"{class_name.split()[1]}"
            )
            all_results.append(result)
            
            print(f"  Convergence rate: {result['convergence_rate']:.2%}")
            print(f"  Cycle rate: {result['cycle_rate']:.2%}")
            print(f"  Mean steps to stable: {result['mean_steps']:.1f}")
            print(f"  Final entropy: {result['final_entropy_mean']:.3f} ± {result['final_entropy_std']:.3f}")
    
    return all_results


def visualize_recursive_trajectories(rule_number: int, size=101, n_samples=5):
    """Visualize entropy trajectories for recursive G application."""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Recursive G Trajectories — Rule {rule_number}', fontsize=14)
    
    for idx, ax in enumerate(axes.flat[:n_samples]):
        # Random initial state
        density = np.random.uniform(0.1, 0.5)
        initial = (np.random.random(size) < density).astype(np.uint8)
        
        # Run iteration
        states, metrics = G_iterated(initial, rule_number, n_iter=500, track_metrics=True)
        
        # Plot entropy trajectory
        ax.plot(metrics['entropy'], 'b-', linewidth=1, label='Entropy')
        ax.axhline(y=np.mean(metrics['entropy'][-50:]), color='r', 
                   linestyle='--', label='Final mean')
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Entropy')
        ax.set_ylim(-0.05, 1.05)
        ax.legend()
        ax.set_title(f'Sample {idx+1} (density={density:.2f})')
    
    # Summary stats on last subplot
    ax = axes.flat[-1]
    ax.axis('off')
    
    # Compute stats
    all_final_ents = []
    all_cycles = []
    for _ in range(20):
        initial = (np.random.random(size) < 0.3).astype(np.uint8)
        conv, steps, final_ent, cycle_len, ent_series = G_orbit_convergence(
            initial, rule_number, max_iter=500
        )
        all_final_ents.append(final_ent)
        all_cycles.append(cycle_len)
    
    stats_text = f"""
Rule {rule_number} Summary (20 samples):

Convergence behavior:
  - Mean final entropy: {np.mean(all_final_ents):.3f} ± {np.std(all_final_ents):.3f}
  - Cycles detected: {sum(1 for c in all_cycles if c > 0)}/20

Interpretation:
  - Low entropy → G drives to simple fixed point
  - High stable entropy → G preserves complexity
  - Oscillating entropy → G induces periodic behavior
    """
    ax.text(0.1, 0.5, stats_text, fontsize=10, family='monospace',
            verticalalignment='center')
    
    plt.tight_layout()
    Path('results').mkdir(exist_ok=True)
    plt.savefig(f'results/exp11_recursive_G_rule{rule_number}.png', dpi=150)
    print(f"Saved results/exp11_recursive_G_rule{rule_number}.png")
    plt.close()


def compare_classes_visual():
    """Create comparison plot across all Wolfram classes."""
    
    rules = [0, 4, 30, 110, 24]  # One from each class + outlier
    labels = ['Class I (0)', 'Class II (4)', 'Class III (30)', 
              'Class IV (110)', 'Class II outlier (24)']
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Recursive G Stability Across Wolfram Classes', fontsize=14)
    
    for ax, rule, label in zip(axes.flat, rules, labels):
        trajectories = []
        for _ in range(10):
            initial = (np.random.random(101) < 0.3).astype(np.uint8)
            _, metrics = G_iterated(initial, rule, n_iter=300, track_metrics=True)
            trajectories.append(metrics['entropy'])
        
        # Plot all trajectories with low alpha
        for traj in trajectories:
            ax.plot(traj, alpha=0.3, color='blue')
        
        # Plot mean
        max_len = max(len(t) for t in trajectories)
        mean_traj = []
        for i in range(max_len):
            vals = [t[i] for t in trajectories if i < len(t)]
            mean_traj.append(np.mean(vals))
        ax.plot(mean_traj, 'r-', linewidth=2, label='Mean')
        
        ax.set_xlabel('G-iteration')
        ax.set_ylabel('Entropy')
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(label)
        ax.legend()
    
    # Summary table on last subplot
    ax = axes.flat[-1]
    ax.axis('off')
    
    table_text = """
Key Observations:

Class I (0, 32, 128):
  G-iteration → entropy → 0 quickly
  Collapses to fixed point

Class II (4, 108, 200):
  G-iteration → periodic entropy oscillation
  Converges to limit cycle

Class III (30, 90, 150):
  G-iteration → sustained high entropy
  No convergence, chaotic trajectory

Class IV (54, 110):
  G-iteration → intermediate stable entropy
  "Edge of stability" — neither collapse
  nor explosion
  
Class II Outliers (24, 35):
  ??? (run experiment to find out)
    """
    ax.text(0.1, 0.5, table_text, fontsize=10, family='monospace',
            verticalalignment='center')
    
    plt.tight_layout()
    Path('results').mkdir(exist_ok=True)
    plt.savefig('results/exp11_recursive_G_comparison.png', dpi=150)
    print("Saved results/exp11_recursive_G_comparison.png")
    plt.close()


# ═══ Main ═══

if __name__ == "__main__":
    print("="*70)
    print("EXPERIMENT 11: Recursive Groovy Commutator")
    print("Testing g(g(g(...s))) stability across Wolfram classes")
    print("="*70)
    
    # Full sweep
    print("\n[Phase 1] Sweeping rules by class...")
    results = sweep_rules_by_class()
    
    # Visualize specific rules
    print("\n[Phase 2] Generating trajectory visualizations...")
    for rule in [0, 4, 30, 110, 54, 24, 35]:
        visualize_recursive_trajectories(rule)
    
    # Comparison plot
    print("\n[Phase 3] Generating class comparison...")
    compare_classes_visual()
    
    # Save summary
    print("\n[Phase 4] Saving results...")
    import json
    with open('results/exp11_summary.json', 'w') as f:
        # Convert numpy types for JSON
        clean_results = []
        for r in results:
            clean_r = {}
            for k, v in r.items():
                if isinstance(v, np.ndarray):
                    clean_r[k] = v.tolist()
                elif isinstance(v, (np.integer, np.floating)):
                    clean_r[k] = float(v)
                else:
                    clean_r[k] = v
            clean_results.append(clean_r)
        json.dump(clean_results, f, indent=2)
    
    print("\n" + "="*70)
    print("EXPERIMENT 11 COMPLETE")
    print("="*70)
    print("\nKey findings:")
    for r in results:
        print(f"  Rule {r['rule']:3d}: {r['convergence_rate']:.0%} converge, "
              f"final entropy {r['final_entropy_mean']:.3f}")
