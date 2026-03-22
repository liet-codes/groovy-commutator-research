"""Elementary Cellular Automata engine using numpy."""

import numpy as np


def rule_lookup(rule_number: int) -> np.ndarray:
    """Return the 8-entry lookup table for an ECA rule (0-255).

    Index i of the table gives the output for the neighborhood pattern
    whose binary representation is i (from 000 to 111).
    """
    return np.array([(rule_number >> i) & 1 for i in range(8)], dtype=np.uint8)


def step(state: np.ndarray, rule_number: int) -> np.ndarray:
    """Advance a 1-D binary state by one ECA step with periodic boundary conditions."""
    table = rule_lookup(rule_number)
    n = len(state)
    # Build 3-bit neighborhood index: left*4 + center*2 + right
    left = np.roll(state, 1)
    right = np.roll(state, -1)
    idx = (left.astype(np.int32) << 2) | (state.astype(np.int32) << 1) | right.astype(np.int32)
    return table[idx]


def evolve(initial: np.ndarray, rule_number: int, timesteps: int) -> np.ndarray:
    """Evolve an ECA for the given number of timesteps.

    Returns an array of shape (timesteps+1, width) where row 0 is the initial state.
    """
    width = len(initial)
    history = np.zeros((timesteps + 1, width), dtype=np.uint8)
    history[0] = initial
    for t in range(timesteps):
        history[t + 1] = step(history[t], rule_number)
    return history


def random_initial(width: int, rng=None) -> np.ndarray:
    """Generate a random binary initial condition."""
    if rng is None:
        rng = np.random.default_rng()
    return rng.integers(0, 2, size=width, dtype=np.uint8)
