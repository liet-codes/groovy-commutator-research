"""Margolus block CA engine for 2-cell binary blocks.

A block rule maps 2-cell blocks (a, b) -> (L(a,b), R(a,b)).
The Margolus scheme alternates even and odd block boundaries:
  - Even phase: blocks at (0,1), (2,3), (4,5), ...
  - Odd phase: blocks at (1,2), (3,4), (5,6), ...

Reference: Brooklyn Rose, "Elementary Interaction Graphs" (2026), §5-6.
"""

import numpy as np


def block_rule_lookup(block_rule_number: int) -> tuple:
    """Decode a block rule number (0-255) into L and R truth tables.
    
    A 2-cell block rule maps (a,b) -> (L(a,b), R(a,b)).
    There are 4 inputs and 4 possible input combinations, so:
      - L is a 2-input Boolean function (4-bit truth table, rules 0-15)
      - R is a 2-input Boolean function (4-bit truth table, rules 0-15)
      - Block rule = L * 16 + R (256 total)
    
    Returns (L_table, R_table) each of shape (4,) indexed by a*2+b.
    """
    R_num = block_rule_number % 16
    L_num = block_rule_number // 16
    
    L_table = np.array([(L_num >> i) & 1 for i in range(4)], dtype=np.uint8)
    R_table = np.array([(R_num >> i) & 1 for i in range(4)], dtype=np.uint8)
    
    return L_table, R_table


def block_step_even(state: np.ndarray, L_table: np.ndarray, R_table: np.ndarray) -> np.ndarray:
    """Apply even-phase block update: blocks at (0,1), (2,3), ..."""
    n = len(state)
    result = state.copy()
    for i in range(0, n - 1, 2):
        a, b = state[i], state[i + 1]
        idx = a * 2 + b
        result[i] = L_table[idx]
        result[i + 1] = R_table[idx]
    return result


def block_step_odd(state: np.ndarray, L_table: np.ndarray, R_table: np.ndarray) -> np.ndarray:
    """Apply odd-phase block update: blocks at (1,2), (3,4), ..."""
    n = len(state)
    result = state.copy()
    for i in range(1, n - 1, 2):
        a, b = state[i], state[i + 1]
        idx = a * 2 + b
        result[i] = L_table[idx]
        result[i + 1] = R_table[idx]
    return result


def block_full_step(state: np.ndarray, L_table: np.ndarray, R_table: np.ndarray) -> np.ndarray:
    """One full Margolus step = even phase + odd phase."""
    intermediate = block_step_even(state, L_table, R_table)
    return block_step_odd(intermediate, L_table, R_table)


def block_evolve(initial: np.ndarray, block_rule_number: int, timesteps: int) -> np.ndarray:
    """Evolve a Margolus block CA for the given number of full steps.
    
    Each full step = even phase + odd phase (2 half-steps).
    Returns array of shape (timesteps+1, width).
    """
    L_table, R_table = block_rule_lookup(block_rule_number)
    width = len(initial)
    history = np.zeros((timesteps + 1, width), dtype=np.uint8)
    history[0] = initial.copy()
    for t in range(timesteps):
        history[t + 1] = block_full_step(history[t], L_table, R_table)
    return history


def block_D(state: np.ndarray, block_rule_number: int) -> np.ndarray:
    """D(S) = S ⊕ φ(S) for block CA. φ = one full Margolus step."""
    L_table, R_table = block_rule_lookup(block_rule_number)
    evolved = block_full_step(state, L_table, R_table)
    return state ^ evolved


def block_G(state: np.ndarray, block_rule_number: int) -> np.ndarray:
    """G(S) = C(D(E(S)), E(D(S))) for block CA.
    
    Path 1: D(E(S)) = E(S) ⊕ E(E(S))
    Path 2: E(D(S)) = φ(D(S))
    """
    L_table, R_table = block_rule_lookup(block_rule_number)
    
    # E(S) = φ(S)
    e_s = block_full_step(state, L_table, R_table)
    
    # Path 1: D(E(S)) = E(S) ⊕ E(E(S))
    e_e_s = block_full_step(e_s, L_table, R_table)
    d_of_e = e_s ^ e_e_s
    
    # Path 2: E(D(S)) = φ(D(S))
    d_s = state ^ e_s  # D(S) = S ⊕ φ(S)
    e_of_d = block_full_step(d_s, L_table, R_table)
    
    # G = C(Path1, Path2) = Path1 ⊕ Path2
    return d_of_e ^ e_of_d


def find_block_rule_30():
    """Find the block rule number for Brooklyn's Block Rule 30.
    
    Block Rule 30: L(x,y) = x, R(x,y) = x ⊕ y.
    
    L truth table (indexed by x*2+y):
      L(0,0)=0, L(0,1)=0, L(1,0)=1, L(1,1)=1 -> L = 0b1100 = 12
    
    R truth table:
      R(0,0)=0, R(0,1)=1, R(1,0)=1, R(1,1)=0 -> R = 0b0110 = 6
    
    Block rule = L*16 + R = 12*16 + 6 = 198
    """
    return 198


def verify_block_rule_30():
    """Verify that block rule 198 gives L(x,y)=x, R(x,y)=x⊕y."""
    L_table, R_table = block_rule_lookup(198)
    
    expected_L = np.array([0, 0, 1, 1], dtype=np.uint8)  # L(x,y) = x
    expected_R = np.array([0, 1, 1, 0], dtype=np.uint8)  # R(x,y) = x⊕y
    
    assert np.array_equal(L_table, expected_L), f"L mismatch: {L_table} vs {expected_L}"
    assert np.array_equal(R_table, expected_R), f"R mismatch: {R_table} vs {expected_R}"
    return True
