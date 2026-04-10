"""
Experiment 12: Glider-Level Groovy Commutator

Hypothesis: Coarse-graining to glider level before applying G will reveal
cleaner convergence for Class IV rules, because we're operating at the
"natural scale" of the system's structure.

Rule 110 glider types (simplified):
- A (right-moving, period 7, speed c/7)
- B (right-moving, period 7, speed c/7, different phase)  
- C (stationary, period 7)
- Various collision products

Approach:
1. Detect gliders in raw CA state
2. Represent as list of (position, type, phase)
3. Define glider-level evolution (collision rules)
4. Apply G at glider level
5. Compare convergence to cell-level G
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import namedtuple

from src.eca import step, evolve
from src.operators import G as cell_level_G
from src.measures import shannon_entropy


# ═══ Glider Definitions for Rule 110 ═══

# Simplified glider patterns (7-cell neighborhoods that characterize gliders)
# These are approximate signatures for detection
GLIDER_SIGNATURES = {
    'A': np.array([0, 1, 1, 1, 0, 0, 0]),  # Right-moving, phase 0
    'A2': np.array([0, 0, 1, 1, 1, 0, 0]),  # Phase variant
    'B': np.array([1, 0, 1, 1, 0, 0, 0]),   # Right-moving, different structure
    'C': np.array([0, 0, 1, 0, 1, 1, 1]),   # Stationary/periodic
}

Glider = namedtuple('Glider', ['position', 'gtype', 'phase'])


# ═══ Glider Detection ═══

def detect_gliders(state: np.ndarray, rule_number: int = 110) -> list:
    """
    Detect gliders in a Rule 110 state.
    
    Returns list of Glider objects. This is a heuristic detector based on
    local pattern matching.
    """
    if rule_number != 110:
        raise ValueError("Glider detection only implemented for Rule 110")
    
    gliders = []
    n = len(state)
    
    # Slide a window and match signatures
    for i in range(n):
        window = np.array([state[(i + j) % n] for j in range(7)])
        
        for gtype, sig in GLIDER_SIGNATURES.items():
            if np.array_equal(window, sig):
                gliders.append(Glider(position=i, gtype=gtype, phase=0))
    
    return gliders


def gliders_to_state(gliders: list, size: int = 101, rule_number: int = 110) -> np.ndarray:
    """
    Convert glider list back to cell state (for visualization/comparison).
    Places glider signatures at their positions.
    """
    state = np.zeros(size, dtype=np.uint8)
    
    for g in gliders:
        sig = GLIDER_SIGNATURES.get(g.gtype, GLIDER_SIGNATURES['A'])
        for j, val in enumerate(sig):
            pos = (g.position + j) % size
            state[pos] = max(state[pos], val)  # Don't overwrite if overlapping
    
    return state


# ═══ Glider-Level Evolution ═══

def evolve_gliders(gliders: list, steps: int = 1, size: int = 101) -> list:
    """
    Evolve gliders forward in time.
    
    This is a simplified model:
    - A gliders move right at c/7 (1 cell per 7 steps, approx 0.14 cells/step)
    - B gliders move right at c/7
    - C gliders stay roughly stationary
    - When gliders get close, they may collide (simplified: they annihilate or transform)
    """
    new_gliders = []
    
    for g in gliders:
        # Simple motion rules
        if g.gtype in ('A', 'A2', 'B'):
            new_pos = (g.position + steps) % size  # Simplified: 1 cell per step
        else:
            new_pos = g.position
        
        new_phase = (g.phase + steps) % 7
        new_gliders.append(Glider(position=new_pos, gtype=g.gtype, phase=new_phase))
    
    # Handle collisions (very simplified: gliders within 3 cells interact)
    if len(new_gliders) > 1:
        # Sort by position
        new_gliders.sort(key=lambda g: g.position)
        
        # Mark collisions
        to_remove = set()
        for i in range(len(new_gliders) - 1):
            for j in range(i + 1, len(new_gliders)):
                dist = min(
                    abs(new_gliders[i].position - new_gliders[j].position),
                    size - abs(new_gliders[i].position - new_gliders[j].position)
                )
                if dist < 5:  # Collision threshold
                    # Simplified: they merge/transform
                    to_remove.add(i)
                    to_remove.add(j)
                    # Could add collision product here
        
        new_gliders = [g for idx, g in enumerate(new_gliders) if idx not in to_remove]
    
    return new_gliders


# ═══ Glider-Level Operators ═══

def D_glider(gliders: list, size: int = 101) -> list:
    """
    Differentiation at glider level: which gliders are about to change?
    
    A glider 'changes' if:
    - It's about to collide with another
    - It's at a phase transition point
    """
    if len(gliders) == 0:
        return []
    
    # Sort by position
    sorted_g = sorted(gliders, key=lambda g: g.position)
    
    changes = []
    for i, g in enumerate(sorted_g):
        # Check for upcoming collision
        next_g = sorted_g[(i + 1) % len(sorted_g)]
        dist = min(
            abs(g.position - next_g.position),
            size - abs(g.position - next_g.position)
        )
        
        # Gliders within 10 cells are 'about to interact'
        if dist < 10:
            changes.append(Glider(position=g.position, gtype=g.gtype, phase=g.phase))
    
    return changes


def I_glider(gliders: list, changes: list) -> list:
    """
    Integration at glider level: apply changes to glider list.
    
    Simplified: remove gliders that are in the change list.
    """
    change_positions = {g.position for g in changes}
    return [g for g in gliders if g.position not in change_positions]


def G_glider(gliders: list, size: int = 101) -> list:
    """
    Groovy Commutator at glider level: [E, D] for gliders.
    
    Path 1: Evolve then differentiate
    Path 2: Differentiate then evolve
    Compare and return difference.
    """
    # Path 1: Evolve then differentiate
    evolved = evolve_gliders(gliders, steps=1, size=size)
    d_of_e = D_glider(evolved, size)
    
    # Path 2: Differentiate then evolve
    d_gliders = D_glider(gliders, size)
    e_of_d = evolve_gliders(d_gliders, steps=1, size=size)
    
    # Commutator = symmetric difference of the two paths
    # Simplified: gliders that appear in one path but not both
    e_set = {(g.position, g.gtype) for g in d_of_e}
    d_set = {(g.position, g.gtype) for g in e_of_d}
    
    sym_diff = e_set.symmetric_difference(d_set)
    
    # Reconstruct glider list from symmetric difference
    result = []
    for g in d_of_e + e_of_d:
        key = (g.position, g.gtype)
        if key in sym_diff:
            result.append(g)
            sym_diff.remove(key)  # Avoid duplicates
    
    return result


# ═══ Measures at Glider Level ═══

def glider_entropy(gliders: list, size: int = 101) -> float:
    """Entropy based on glider density and type distribution."""
    if len(gliders) == 0:
        return 0.0
    
    # Position entropy (how spread out?)
    positions = [g.position / size for g in gliders]
    pos_entropy = shannon_entropy(np.array(positions) > 0.5)  # Simplified
    
    # Type entropy
    type_counts = {}
    for g in gliders:
        type_counts[g.gtype] = type_counts.get(g.gtype, 0) + 1
    
    total = len(gliders)
    probs = np.array([c / total for c in type_counts.values()])
    type_entropy = -np.sum(probs * np.log2(probs + 1e-10))
    
    # Combined
    return (len(gliders) / (size / 7)) * 0.5 + type_entropy / 2  # Normalize roughly


def glider_count(gliders: list) -> int:
    """Simple count of gliders."""
    return len(gliders)


# ═══ Experiments ═══

def compare_cell_vs_glider_G(rule_number: int = 110, n_samples: int = 10, 
                              n_iter: int = 100, size: int = 101):
    """
    Compare G-iteration at cell level vs glider level for Rule 110.
    """
    if rule_number != 110:
        print(f"Warning: Glider detection only valid for Rule 110, not {rule_number}")
        return
    
    results = {
        'cell_level': [],
        'glider_level': [],
    }
    
    for sample in range(n_samples):
        # Random initial state with some structure
        initial = np.zeros(size, dtype=np.uint8)
        # Seed with some glider-like patterns
        for _ in range(5):
            pos = np.random.randint(0, size - 7)
            gtype = np.random.choice(['A', 'B', 'C'])
            sig = GLIDER_SIGNATURES.get(gtype, GLIDER_SIGNATURES['A'])
            initial[pos:pos+7] = sig
        
        # Cell-level trajectory
        cell_entropies = []
        s = initial.copy()
        for _ in range(n_iter):
            cell_entropies.append(shannon_entropy(s))
            s = cell_level_G(s, rule_number)
        results['cell_level'].append(cell_entropies)
        
        # Glider-level trajectory
        glider_entropies = []
        gliders = detect_gliders(initial, rule_number)
        for _ in range(n_iter):
            glider_entropies.append(glider_entropy(gliders, size))
            gliders = G_glider(gliders, size)
        results['glider_level'].append(glider_entropies)
    
    return results


def visualize_comparison():
    """Visualize cell-level vs glider-level G convergence."""
    
    print("Running cell-level vs glider-level G comparison...")
    results = compare_cell_vs_glider_G(n_samples=15, n_iter=150)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Rule 110: Cell-Level vs Glider-Level G Convergence', fontsize=14)
    
    # Cell level
    ax = axes[0]
    for traj in results['cell_level']:
        ax.plot(traj, alpha=0.3, color='blue', linewidth=0.8)
    mean_cell = np.mean(results['cell_level'], axis=0)
    ax.plot(mean_cell, color='darkblue', linewidth=2, label='Mean')
    ax.axhline(y=mean_cell[-1], color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('G-iteration')
    ax.set_ylabel('Entropy')
    ax.set_ylim(-0.05, 1.05)
    ax.set_title('Cell-Level G')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Glider level
    ax = axes[1]
    for traj in results['glider_level']:
        ax.plot(traj, alpha=0.3, color='green', linewidth=0.8)
    mean_glider = np.mean(results['glider_level'], axis=0)
    ax.plot(mean_glider, color='darkgreen', linewidth=2, label='Mean')
    ax.axhline(y=mean_glider[-1], color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('G-iteration')
    ax.set_ylabel('Glider Entropy')
    ax.set_title('Glider-Level G')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Add text summary
    fig.text(0.5, 0.02, 
             f'Cell-level final entropy: {mean_cell[-1]:.3f} | '
             f'Glider-level final entropy: {mean_glider[-1]:.3f}',
             ha='center', fontsize=11)
    
    plt.tight_layout()
    plt.savefig('/tmp/glider_level_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved /tmp/glider_level_comparison.png")
    plt.close()


def test_glider_stability():
    """Test if glider-level G converges faster/more cleanly."""
    
    print("\nTesting glider-level G stability...")
    
    size = 101
    n_samples = 20
    
    cell_convergence_steps = []
    glider_convergence_steps = []
    
    for sample in range(n_samples):
        # Create initial state with gliders
        initial = np.zeros(size, dtype=np.uint8)
        for _ in range(5):
            pos = np.random.randint(0, size - 7)
            gtype = np.random.choice(['A', 'B', 'C'])
            sig = GLIDER_SIGNATURES.get(gtype, GLIDER_SIGNATURES['A'])
            initial[pos:pos+7] = sig
        
        # Cell-level convergence
        s = initial.copy()
        prev_ent = shannon_entropy(s)
        cell_stable_for = 0
        cell_converged_at = None
        
        for step in range(200):
            s = cell_level_G(s, 110)
            ent = shannon_entropy(s)
            if abs(ent - prev_ent) < 0.01:
                cell_stable_for += 1
                if cell_stable_for >= 20 and cell_converged_at is None:
                    cell_converged_at = step
            else:
                cell_stable_for = 0
            prev_ent = ent
        
        if cell_converged_at:
            cell_convergence_steps.append(cell_converged_at)
        
        # Glider-level convergence
        gliders = detect_gliders(initial, 110)
        prev_g_ent = glider_entropy(gliders, size)
        glider_stable_for = 0
        glider_converged_at = None
        
        for step in range(200):
            gliders = G_glider(gliders, size)
            g_ent = glider_entropy(gliders, size)
            if abs(g_ent - prev_g_ent) < 0.05:
                glider_stable_for += 1
                if glider_stable_for >= 20 and glider_converged_at is None:
                    glider_converged_at = step
            else:
                glider_stable_for = 0
            prev_g_ent = g_ent
        
        if glider_converged_at:
            glider_convergence_steps.append(glider_converged_at)
    
    print(f"\nCell-level convergence: {len(cell_convergence_steps)}/{n_samples} samples")
    if cell_convergence_steps:
        print(f"  Mean steps to converge: {np.mean(cell_convergence_steps):.1f} ± {np.std(cell_convergence_steps):.1f}")
    
    print(f"\nGlider-level convergence: {len(glider_convergence_steps)}/{n_samples} samples")
    if glider_convergence_steps:
        print(f"  Mean steps to converge: {np.mean(glider_convergence_steps):.1f} ± {np.std(glider_convergence_steps):.1f}")
    
    if cell_convergence_steps and glider_convergence_steps:
        ratio = np.mean(cell_convergence_steps) / np.mean(glider_convergence_steps)
        print(f"\nGlider-level G converges {ratio:.2f}x faster (on average)")


# ═══ Main ═══

if __name__ == "__main__":
    print("="*70)
    print("EXPERIMENT 12: Glider-Level Groovy Commutator")
    print("="*70)
    
    # Visual comparison
    visualize_comparison()
    
    # Stability test
    test_glider_stability()
    
    print("\n" + "="*70)
    print("EXPERIMENT 12 COMPLETE")
    print("="*70)
