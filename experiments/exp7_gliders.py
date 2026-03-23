"""
Experiment 7: Glider Detection via Groovy Commutator

Focus on Rule 110 (proven universal).
- Run N=201, T=500
- Compute G at each step
- Track local maxima of |G| that persist across timesteps
- A "glider" should appear as a localized, persistent, moving peak in G
- Generate spacetime diagram of G overlaid with detected glider tracks

Output: results/glider_detection.png, results/glider_stats.csv
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from src.eca import step, random_initial
from src.operators import G as groovy_commutator
from src.measures import shannon_entropy


def find_local_maxima(arr, min_width=2):
    """Find positions of local maxima in a binary/integer array.

    Returns list of (center_position, width) tuples for contiguous nonzero regions.
    """
    regions = []
    n = len(arr)
    i = 0
    while i < n:
        if arr[i] > 0:
            start = i
            while i < n and arr[i] > 0:
                i += 1
            width = i - start
            center = (start + i - 1) / 2.0
            regions.append((center, width))
        else:
            i += 1
    return regions


def track_gliders(g_history, max_speed=3, min_persistence=5):
    """Track persistent moving peaks in G spacetime.

    A glider is a localized peak that:
    - Persists for at least min_persistence timesteps
    - Moves with approximately constant velocity (|speed| <= max_speed cells/step)

    Returns list of glider tracks: [(t_start, positions_list, speed)]
    """
    T, N = g_history.shape
    active_tracks = []   # [(t_start, [positions], last_pos)]
    finished_tracks = []

    for t in range(T):
        peaks = find_local_maxima(g_history[t])
        peak_positions = [p[0] for p in peaks]

        # Try to extend existing tracks
        used_peaks = set()
        new_active = []

        for t_start, positions, last_pos in active_tracks:
            best_match = None
            best_dist = max_speed + 1

            for idx, pp in enumerate(peak_positions):
                if idx in used_peaks:
                    continue
                # Account for periodic boundaries
                dist = min(abs(pp - last_pos), N - abs(pp - last_pos))
                if dist < best_dist:
                    best_dist = dist
                    best_match = idx

            if best_match is not None and best_dist <= max_speed:
                used_peaks.add(best_match)
                positions.append(peak_positions[best_match])
                new_active.append((t_start, positions, peak_positions[best_match]))
            else:
                # Track ended
                if len(positions) >= min_persistence:
                    finished_tracks.append((t_start, positions))

        active_tracks = new_active

        # Start new tracks from unmatched peaks
        for idx, pp in enumerate(peak_positions):
            if idx not in used_peaks:
                active_tracks.append((t, [pp], pp))

    # Finish remaining active tracks
    for t_start, positions, _ in active_tracks:
        if len(positions) >= min_persistence:
            finished_tracks.append((t_start, positions))

    # Compute speeds
    glider_tracks = []
    for t_start, positions in finished_tracks:
        if len(positions) >= 2:
            # Compute velocity accounting for periodic boundary
            velocities = []
            for i in range(1, len(positions)):
                dp = positions[i] - positions[i - 1]
                # Unwrap periodic boundary
                if dp > N / 2:
                    dp -= N
                elif dp < -N / 2:
                    dp += N
                velocities.append(dp)
            speed = np.mean(velocities) if velocities else 0.0
        else:
            speed = 0.0
        glider_tracks.append((t_start, positions, speed))

    return glider_tracks


def run():
    N = 201
    T = 500
    SEED = 42
    RULE = 110

    print(f"Experiment 7: Glider Detection via Groovy Commutator")
    print(f"Rule {RULE}, N={N}, T={T}")

    rng = np.random.default_rng(SEED)
    init = random_initial(N, rng)

    # Evolve and collect G at each step
    state_history = np.zeros((T + 1, N), dtype=np.uint8)
    g_history = np.zeros((T, N), dtype=np.uint8)
    state_history[0] = init

    for t in range(T):
        g = groovy_commutator(state_history[t], RULE)
        g_history[t] = g
        state_history[t + 1] = step(state_history[t], RULE)

    # Also run for Rules 54 and 30 for comparison
    comparison_rules = [54, 30]
    comparison_g = {}
    comparison_states = {}
    for rule in comparison_rules:
        state = random_initial(N, np.random.default_rng(SEED))
        sh = np.zeros((T + 1, N), dtype=np.uint8)
        gh = np.zeros((T, N), dtype=np.uint8)
        sh[0] = state
        for t in range(T):
            gh[t] = groovy_commutator(sh[t], rule)
            sh[t + 1] = step(sh[t], rule)
        comparison_g[rule] = gh
        comparison_states[rule] = sh

    # Detect gliders in Rule 110
    glider_tracks = track_gliders(g_history, max_speed=3, min_persistence=8)
    print(f"\nDetected {len(glider_tracks)} glider tracks in Rule {RULE}")

    # Classify gliders by speed
    if glider_tracks:
        speeds = [abs(t[2]) for t in glider_tracks]
        durations = [len(t[1]) for t in glider_tracks]
        print(f"  Speed range: {min(speeds):.2f} — {max(speeds):.2f} cells/step")
        print(f"  Duration range: {min(durations)} — {max(durations)} steps")

        # Group by approximate speed
        speed_bins = {}
        for t_start, positions, speed in glider_tracks:
            s_round = round(abs(speed), 1)
            if s_round not in speed_bins:
                speed_bins[s_round] = 0
            speed_bins[s_round] += 1
        print("  Speed distribution:")
        for s, count in sorted(speed_bins.items()):
            print(f"    |v|≈{s:.1f}: {count} tracks")

    # Also detect for comparison rules
    for rule in comparison_rules:
        gt = track_gliders(comparison_g[rule], max_speed=3, min_persistence=8)
        print(f"  Rule {rule}: {len(gt)} tracks detected")

    # Write stats CSV
    os.makedirs('results', exist_ok=True)
    stats_rows = []
    for t_start, positions, speed in glider_tracks:
        stats_rows.append({
            'rule': RULE,
            't_start': t_start,
            'duration': len(positions),
            'speed': round(speed, 4),
            'start_position': round(positions[0], 1),
            'end_position': round(positions[-1], 1),
        })

    if stats_rows:
        with open('results/glider_stats.csv', 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=stats_rows[0].keys())
            writer.writeheader()
            writer.writerows(stats_rows)
        print(f"\nSaved results/glider_stats.csv ({len(stats_rows)} tracks)")
    else:
        with open('results/glider_stats.csv', 'w', newline='') as f:
            f.write("rule,t_start,duration,speed,start_position,end_position\n")
        print("\nSaved results/glider_stats.csv (no tracks found)")

    # Generate plots
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    # Row 1: Spacetime diagrams of the CA state
    for ax, (rule, label) in zip(axes[0], [(RULE, f'Rule {RULE} (Class IV)'),
                                            (54, 'Rule 54 (Class IV)'),
                                            (30, 'Rule 30 (Class III)')]):
        if rule == RULE:
            sh = state_history[:-1]
        else:
            sh = comparison_states[rule][:-1]
        ax.imshow(sh[:200], cmap='binary', aspect='auto', interpolation='nearest')
        ax.set_title(f'{label} — CA State', fontsize=10)
        ax.set_xlabel('Position')
        ax.set_ylabel('Time')

    # Row 2: G spacetime with glider tracks overlaid
    for ax, (rule, label) in zip(axes[1], [(RULE, f'Rule {RULE} (Class IV)'),
                                            (54, 'Rule 54 (Class IV)'),
                                            (30, 'Rule 30 (Class III)')]):
        if rule == RULE:
            gh = g_history
            tracks = glider_tracks
        else:
            gh = comparison_g[rule]
            tracks = track_gliders(gh, max_speed=3, min_persistence=8)

        ax.imshow(gh[:200], cmap='hot', aspect='auto', interpolation='nearest')

        # Overlay glider tracks
        colors = plt.cm.Set1(np.linspace(0, 1, max(len(tracks), 1)))
        for i, (t_start, positions, speed) in enumerate(tracks):
            # Only show tracks within displayed time range
            vis_t = []
            vis_p = []
            for j, pos in enumerate(positions):
                t = t_start + j
                if t < 200:
                    vis_t.append(t)
                    vis_p.append(pos)
            if vis_t:
                ax.plot(vis_p, vis_t, '-', color=colors[i % len(colors)],
                        linewidth=1.5, alpha=0.8)

        n_tracks = len(tracks)
        ax.set_title(f'{label} — G with tracks ({n_tracks} detected)', fontsize=10)
        ax.set_xlabel('Position')
        ax.set_ylabel('Time')

    plt.suptitle('Glider Detection via Groovy Commutator', fontsize=14)
    plt.tight_layout()
    plt.savefig('results/glider_detection.png', dpi=150)
    print("Saved results/glider_detection.png")
    plt.close()


if __name__ == '__main__':
    run()
