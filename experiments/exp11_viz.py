"""Quick visualization of recursive G results."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from src.eca import step, evolve
from src.operators import G, D, E
from src.measures import shannon_entropy


def G_iterated(initial_state, rule_number, n_iter):
    """Apply G iteratively."""
    s = initial_state.copy()
    entropies = [shannon_entropy(s)]
    for _ in range(n_iter):
        s = G(s, rule_number)
        entropies.append(shannon_entropy(s))
    return entropies


def plot_comparison():
    """Plot entropy trajectories for key rules."""
    
    rules = [
        (0, 'Class I (0)'),
        (4, 'Class II (4)'), 
        (30, 'Class III (30)'),
        (110, 'Class IV (110)'),
        (54, 'Class IV (54)'),
        (24, 'Class II outlier (24)'),
    ]
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Recursive G Stability: Entropy Trajectories by Wolfram Class', fontsize=14)
    
    colors = {'Class I': '#1f77b4', 'Class II': '#ff7f0e', 'Class III': '#d62728', 
              'Class IV': '#2ca02c', 'Class II outlier': '#9467bd'}
    
    for idx, (rule, label) in enumerate(rules):
        ax = axes.flat[idx]
        class_name = label.split()[0] + ' ' + label.split()[1].strip('()')
        
        # Run multiple samples
        all_trajs = []
        for _ in range(15):
            density = np.random.uniform(0.1, 0.5)
            initial = (np.random.random(101) < density).astype(np.uint8)
            traj = G_iterated(initial, rule, 300)
            all_trajs.append(traj)
        
        # Plot individual trajectories
        for traj in all_trajs:
            ax.plot(traj, alpha=0.25, color=colors.get(class_name, 'gray'), linewidth=0.8)
        
        # Plot mean
        max_len = max(len(t) for t in all_trajs)
        mean_traj = []
        for i in range(max_len):
            vals = [t[i] for t in all_trajs if i < len(t)]
            mean_traj.append(np.mean(vals))
        ax.plot(mean_traj, color='black', linewidth=2.5, label='Mean')
        
        # Styling
        ax.set_xlabel('G-iteration', fontsize=10)
        ax.set_ylabel('Entropy', fontsize=10)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlim(0, 300)
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add final entropy annotation
        final_ent = mean_traj[-1] if mean_traj else 0
        ax.axhline(y=final_ent, color='red', linestyle='--', alpha=0.5)
        ax.text(250, final_ent + 0.05, f'final: {final_ent:.3f}', fontsize=9, color='red')
    
    plt.tight_layout()
    plt.savefig('/tmp/recursive_g_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved /tmp/recursive_g_comparison.png")
    plt.close()


if __name__ == "__main__":
    plot_comparison()
