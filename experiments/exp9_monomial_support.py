"""Experiment 9: Monomial Support as Dynamical Classifier

Tests Ben Sachs' conjecture that degree-profile equivalence predicts
dynamical behavior (specifically Wolfram class and G entropy) as well as
or better than orbit equivalence.

Protocol:
1. Compute ANF and degree profile for all 256 ECA rules
2. Group the 88 equivalence classes by degree profile
3. Compute within-group vs between-group variance of G entropy
4. Test perturbation sensitivity by monomial degree
5. Correlate degree profile features with Wolfram class

Reference: Ben Sachs, "Monomial Support as a Dynamical Classifier" (2026).
"""

import csv
import os
import sys
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.anf import compute_anf, compute_all_anfs, degree_profile_str, anf_polynomial_str, verify_anf
from src.eca import random_initial, evolve, step
from src.operators import G, D, E, phi
from src.measures import shannon_entropy


# ═══ Verification ═══

def verify_all_anfs():
    """Verify ANF computation for all 256 rules."""
    print("Verifying ANF computation for all 256 rules...")
    failures = []
    for r in range(256):
        if not verify_anf(r):
            failures.append(r)
    if failures:
        print(f"  FAILED: rules {failures}")
        return False
    print("  All 256 rules verified ✓")
    return True


# ═══ Load equivalence class data ═══

def load_equivalence_data(filepath):
    """Load equivalence class data from CSV."""
    classes = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            classes.append({
                'representative': int(row['representative']),
                'equivalents': [int(x) for x in row['equivalents'].split('|')],
                'wolfram_class': row['wolfram_class'],
                'g_entropy_mean': float(row['g_entropy_mean']),
                'perturbation_weight': float(row['perturbation_weight']),
                'affine_base': row['affine_base'],
            })
    return classes


# ═══ Analysis ═══

def analyze_degree_profiles(all_anfs, equiv_classes):
    """Group equivalence classes by degree profile, analyze G entropy variance."""
    
    # Map rule -> ANF info
    rule_to_anf = {a['rule']: a for a in all_anfs}
    
    # Add ANF info to equivalence classes
    for cls in equiv_classes:
        rep = cls['representative']
        anf_info = rule_to_anf[rep]
        cls['degree_profile'] = anf_info['degree_profile']
        cls['algebraic_degree'] = anf_info['algebraic_degree']
        cls['n_monomials'] = anf_info['n_monomials']
        cls['support'] = anf_info['support']
        cls['anf_polynomial'] = anf_polynomial_str(anf_info['anf_coeffs'])
    
    # Group by degree profile
    profile_groups = defaultdict(list)
    for cls in equiv_classes:
        profile_groups[cls['degree_profile']].append(cls)
    
    return profile_groups


def compute_variance_ratio(profile_groups):
    """Compute within-group vs between-group variance ratio for G entropy.
    
    If degree profile is a good predictor, within-group variance should be
    small relative to between-group variance (low ratio = good predictor).
    """
    # Filter to groups with nonzero G entropy (skip trivial rules)
    all_entropies = []
    group_means = []
    within_vars = []
    
    for profile, classes in profile_groups.items():
        entropies = [c['g_entropy_mean'] for c in classes]
        all_entropies.extend(entropies)
        if len(entropies) > 1:
            group_means.append(np.mean(entropies))
            within_vars.append(np.var(entropies))
    
    total_variance = np.var(all_entropies) if all_entropies else 0
    mean_within_variance = np.mean(within_vars) if within_vars else 0
    between_variance = np.var(group_means) if group_means else 0
    
    return {
        'total_variance': total_variance,
        'mean_within_variance': mean_within_variance,
        'between_variance': between_variance,
        'ratio': mean_within_variance / total_variance if total_variance > 0 else float('inf'),
        'n_groups': len(profile_groups),
        'n_nontrivial_groups': len([g for g in profile_groups.values() if len(g) > 1]),
    }


def wolfram_class_purity(profile_groups):
    """For each degree profile group, measure Wolfram class purity.
    
    Purity = fraction of classes in the group that share the most common
    Wolfram class. High purity = degree profile predicts Wolfram class.
    """
    results = []
    for profile, classes in sorted(profile_groups.items()):
        class_counts = defaultdict(int)
        for c in classes:
            class_counts[c['wolfram_class']] += 1
        total = len(classes)
        majority_class = max(class_counts, key=class_counts.get)
        purity = class_counts[majority_class] / total
        results.append({
            'profile': profile,
            'n_classes': total,
            'majority_class': majority_class,
            'purity': purity,
            'class_distribution': dict(class_counts),
            'g_entropy_range': (
                min(c['g_entropy_mean'] for c in classes),
                max(c['g_entropy_mean'] for c in classes),
            ),
        })
    return results


def perturbation_by_degree(all_anfs):
    """Analyze single-bit perturbation effects by monomial degree.
    
    For each rule, flip each output bit (changing support by one monomial),
    measure the change in G entropy. Group results by degree of the changed
    monomial.
    """
    rng = np.random.default_rng(42)
    N = 101
    T = 200
    n_samples = 20
    
    results_by_degree = defaultdict(list)
    
    rule_to_anf = {a['rule']: a for a in all_anfs}
    
    # Sample a diverse set of rules (all Class III and IV, plus some Class II)
    test_rules = [18, 22, 30, 45, 54, 60, 90, 105, 106, 110, 122, 126, 130, 146, 150]
    
    print(f"\nPerturbation analysis: {len(test_rules)} rules × 8 flips × {n_samples} samples...")
    
    for rule in test_rules:
        # Baseline G entropy
        base_entropies = []
        for _ in range(n_samples):
            state = rng.integers(0, 2, size=N, dtype=np.uint8)
            g_ents = []
            for t in range(T):
                g_state = G(state, rule)
                g_ents.append(shannon_entropy(g_state))
                state = phi(state, rule)
            base_entropies.append(np.mean(g_ents[50:]))  # skip transient
        
        base_g = np.mean(base_entropies)
        
        # Now flip each output bit (= change one monomial in ANF)
        base_anf = rule_to_anf[rule]
        
        for bit_pos in range(8):
            perturbed_rule = rule ^ (1 << bit_pos)
            perturbed_anf = rule_to_anf[perturbed_rule]
            
            # Which monomial changed?
            changed_monomials = base_anf['support'].symmetric_difference(perturbed_anf['support'])
            # In GF(2), flipping one truth table bit changes exactly one ANF coefficient
            # But due to the Möbius transform, flipping truth table bit k can affect
            # multiple ANF coefficients. Let's compute it directly.
            degree_changed = perturbed_anf['algebraic_degree'] - base_anf['algebraic_degree']
            
            # Measure perturbed G entropy
            perturbed_entropies = []
            for _ in range(n_samples):
                state = rng.integers(0, 2, size=N, dtype=np.uint8)
                g_ents = []
                for t in range(T):
                    g_state = G(state, perturbed_rule)
                    g_ents.append(shannon_entropy(g_state))
                    state = phi(state, perturbed_rule)
                perturbed_entropies.append(np.mean(g_ents[50:]))
            
            perturbed_g = np.mean(perturbed_entropies)
            delta_g = abs(perturbed_g - base_g)
            
            # Classify by what changed in the support
            support_diff = base_anf['support'].symmetric_difference(perturbed_anf['support'])
            if support_diff:
                max_degree_changed = max(
                    sum(1 for b in range(3) if idx & (1 << b)) for idx in support_diff
                )
            else:
                max_degree_changed = 0
            
            results_by_degree[max_degree_changed].append({
                'rule': rule,
                'bit_flipped': bit_pos,
                'perturbed_rule': perturbed_rule,
                'base_g': base_g,
                'perturbed_g': perturbed_g,
                'delta_g': delta_g,
                'n_support_changes': len(support_diff),
            })
    
    return results_by_degree


def main():
    print("=" * 70)
    print("EXPERIMENT 9: Monomial Support as Dynamical Classifier")
    print("Testing Ben Sachs' conjecture (March 2026)")
    print("=" * 70)
    
    # Step 0: Verify
    if not verify_all_anfs():
        print("ANF verification failed! Aborting.")
        return
    
    # Step 1: Compute all ANFs
    print("\n--- Step 1: Computing ANF for all 256 rules ---")
    all_anfs = compute_all_anfs()
    
    # Print the ANF table
    print(f"\n{'Rule':>5} {'Polynomial':>35} {'Deg Profile':>12} {'Alg Deg':>8}")
    print("-" * 65)
    for a in all_anfs[:20]:  # First 20
        print(f"{a['rule']:>5} {anf_polynomial_str(a['anf_coeffs']):>35} "
              f"{degree_profile_str(a['degree_profile']):>12} {a['algebraic_degree']:>8}")
    print(f"  ... ({len(all_anfs)} total)")
    
    # Count distinct degree profiles
    profiles = set(a['degree_profile'] for a in all_anfs)
    print(f"\nDistinct degree profiles across all 256 rules: {len(profiles)}")
    for p in sorted(profiles):
        count = sum(1 for a in all_anfs if a['degree_profile'] == p)
        print(f"  {degree_profile_str(p):>12}: {count} rules")
    
    # Step 2: Load equivalence class data and group by degree profile
    print("\n--- Step 2: Grouping equivalence classes by degree profile ---")
    equiv_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              'results', 'equivalence_classes.csv')
    equiv_classes = load_equivalence_data(equiv_path)
    profile_groups = analyze_degree_profiles(all_anfs, equiv_classes)
    
    print(f"\n{len(profile_groups)} distinct degree profiles across 88 equivalence classes:")
    for profile, classes in sorted(profile_groups.items()):
        wolfram_classes = [c['wolfram_class'] for c in classes]
        g_entropies = [c['g_entropy_mean'] for c in classes]
        class_counts = defaultdict(int)
        for wc in wolfram_classes:
            class_counts[wc] += 1
        print(f"  {degree_profile_str(profile):>12}: {len(classes):>2} equiv classes | "
              f"G entropy: [{min(g_entropies):.3f}, {max(g_entropies):.3f}] | "
              f"Wolfram: {dict(class_counts)}")
    
    # Step 3: Variance analysis
    print("\n--- Step 3: Variance analysis ---")
    vr = compute_variance_ratio(profile_groups)
    print(f"  Total variance of G entropy: {vr['total_variance']:.6f}")
    print(f"  Mean within-group variance:  {vr['mean_within_variance']:.6f}")
    print(f"  Between-group variance:      {vr['between_variance']:.6f}")
    print(f"  Within/Total ratio:          {vr['ratio']:.4f}")
    print(f"  (Lower ratio = degree profile explains more variance)")
    print(f"  Number of degree profile groups: {vr['n_groups']}")
    
    # Step 4: Wolfram class purity
    print("\n--- Step 4: Wolfram class purity by degree profile ---")
    purity_results = wolfram_class_purity(profile_groups)
    for pr in purity_results:
        print(f"  {degree_profile_str(pr['profile']):>12}: {pr['n_classes']:>2} classes, "
              f"purity={pr['purity']:.2f} (majority: {pr['majority_class']}) "
              f"G range=[{pr['g_entropy_range'][0]:.3f}, {pr['g_entropy_range'][1]:.3f}] "
              f"{pr['class_distribution']}")
    
    mean_purity = np.mean([pr['purity'] for pr in purity_results])
    weighted_purity = sum(pr['purity'] * pr['n_classes'] for pr in purity_results) / 88
    print(f"\n  Mean purity: {mean_purity:.3f}")
    print(f"  Weighted purity: {weighted_purity:.3f}")
    print(f"  (1.0 = degree profile perfectly predicts Wolfram class)")
    
    # Step 5: ANF details for Class IV rules
    print("\n--- Step 5: Class IV rule ANF details ---")
    class_iv = [c for c in equiv_classes if c['wolfram_class'] == 'IV']
    for c in class_iv:
        print(f"  Rule {c['representative']}: {c['anf_polynomial']}")
        print(f"    Degree profile: {degree_profile_str(c['degree_profile'])}, "
              f"Algebraic degree: {c['algebraic_degree']}, "
              f"G entropy: {c['g_entropy_mean']:.4f}")
    
    # Step 6: Perturbation analysis by degree
    print("\n--- Step 6: Perturbation sensitivity by monomial degree ---")
    perturb_results = perturbation_by_degree(all_anfs)
    
    print(f"\n  Perturbation results by max degree of changed monomial:")
    for degree in sorted(perturb_results.keys()):
        deltas = [r['delta_g'] for r in perturb_results[degree]]
        print(f"  Degree {degree}: n={len(deltas):>3}, "
              f"mean |ΔG| = {np.mean(deltas):.4f} ± {np.std(deltas):.4f}, "
              f"max |ΔG| = {np.max(deltas):.4f}")
    
    # Step 7: Save results
    print("\n--- Step 7: Saving results ---")
    results_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results')
    
    # Save ANF table
    anf_path = os.path.join(results_dir, 'anf_table.csv')
    with open(anf_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rule', 'polynomial', 'degree_profile', 'algebraic_degree',
                         'n_monomials', 'support_indices'])
        for a in all_anfs:
            writer.writerow([
                a['rule'],
                anf_polynomial_str(a['anf_coeffs']),
                degree_profile_str(a['degree_profile']),
                a['algebraic_degree'],
                a['n_monomials'],
                '|'.join(str(x) for x in sorted(a['support'])),
            ])
    print(f"  Saved ANF table: {anf_path}")
    
    # Save degree profile analysis
    dp_path = os.path.join(results_dir, 'degree_profile_analysis.csv')
    with open(dp_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['degree_profile', 'n_equiv_classes', 'wolfram_purity',
                         'majority_class', 'g_entropy_min', 'g_entropy_max',
                         'class_distribution'])
        for pr in purity_results:
            writer.writerow([
                degree_profile_str(pr['profile']),
                pr['n_classes'],
                f"{pr['purity']:.3f}",
                pr['majority_class'],
                f"{pr['g_entropy_range'][0]:.4f}",
                f"{pr['g_entropy_range'][1]:.4f}",
                str(pr['class_distribution']),
            ])
    print(f"  Saved degree profile analysis: {dp_path}")
    
    # Save perturbation analysis
    perturb_path = os.path.join(results_dir, 'perturbation_by_degree.csv')
    with open(perturb_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['max_degree_changed', 'n_perturbations', 'mean_delta_g',
                         'std_delta_g', 'max_delta_g'])
        for degree in sorted(perturb_results.keys()):
            deltas = [r['delta_g'] for r in perturb_results[degree]]
            writer.writerow([
                degree,
                len(deltas),
                f"{np.mean(deltas):.6f}",
                f"{np.std(deltas):.6f}",
                f"{np.max(deltas):.6f}",
            ])
    print(f"  Saved perturbation analysis: {perturb_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Distinct degree profiles (all 256 rules): {len(profiles)}")
    print(f"  Distinct degree profiles (88 equiv classes): {len(profile_groups)}")
    print(f"  Variance ratio (within/total): {vr['ratio']:.4f}")
    print(f"  Weighted Wolfram class purity: {weighted_purity:.3f}")
    print(f"  Class IV degree profiles: {[degree_profile_str(c['degree_profile']) for c in class_iv]}")
    print()
    if vr['ratio'] < 0.5:
        print("  → Degree profile explains >50% of G entropy variance")
        print("    Ben's conjecture SUPPORTED at the ECA level")
    else:
        print("  → Degree profile explains <50% of G entropy variance")
        print("    Ben's conjecture NOT SUPPORTED at the ECA level")
        print("    (May still hold for larger rule spaces with more degree profiles)")


if __name__ == '__main__':
    main()
