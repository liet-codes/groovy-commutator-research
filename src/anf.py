"""Algebraic Normal Form (ANF) for ECA rules.

Computes the unique ANF representation of each ECA rule as a polynomial
over GF(2). For n=3 inputs (L, C, R), the 8 possible monomials are:

  degree 0: 1 (constant)
  degree 1: L, C, R
  degree 2: LC, LR, CR
  degree 3: LCR

The ANF coefficients are computed from the truth table via the Möbius
(butterfly) transform over GF(2).

Reference: Ben Sachs, "Monomial Support as a Dynamical Classifier" (2026).
"""

import numpy as np
from src.eca import rule_lookup


# Monomial labels indexed by their binary mask over (L, C, R)
# Bit 0 = L (leftmost input), Bit 1 = C (center), Bit 2 = R (rightmost)
MONOMIAL_NAMES = {
    0b000: "1",      # constant
    0b001: "L",      # degree 1
    0b010: "C",      # degree 1
    0b100: "R",      # degree 1
    0b011: "LC",     # degree 2
    0b101: "LR",     # degree 2
    0b110: "CR",     # degree 2
    0b111: "LCR",    # degree 3
}

MONOMIAL_DEGREES = {
    0b000: 0,
    0b001: 1,
    0b010: 1,
    0b100: 1,
    0b011: 2,
    0b101: 2,
    0b110: 2,
    0b111: 3,
}


def mobius_transform(truth_table: np.ndarray) -> np.ndarray:
    """Compute ANF coefficients from truth table via Möbius butterfly transform.
    
    For n variables, the truth table has 2^n entries. The transform computes
    the unique polynomial over GF(2) that matches the truth table.
    
    This is an in-place butterfly over GF(2): for each variable i,
    XOR entries whose index differs in bit i.
    """
    n = len(truth_table)
    anf = truth_table.copy().astype(np.uint8)
    bit = 1
    while bit < n:
        for i in range(n):
            if i & bit:
                anf[i] ^= anf[i ^ bit]
        bit <<= 1
    return anf


def compute_anf(rule_number: int) -> dict:
    """Compute the full ANF decomposition of an ECA rule.
    
    Returns dict with:
      - 'rule': rule number
      - 'truth_table': 8-entry lookup table
      - 'anf_coeffs': 8-entry ANF coefficient vector (over GF(2))
      - 'support': set of monomial indices with nonzero coefficient
      - 'support_names': set of monomial names with nonzero coefficient
      - 'degree_profile': tuple (d0, d1, d2, d3) counting monomials at each degree
      - 'algebraic_degree': maximum degree of any monomial in the support
      - 'n_monomials': total number of nonzero monomials
    """
    table = rule_lookup(rule_number)
    # Note: rule_lookup indexes by neighborhood pattern where
    # index i = L*4 + C*2 + R. But for Möbius transform, we need
    # the truth table indexed by (x0, x1, x2) = (L, C, R) where
    # index j = L*1 + C*2 + R*4.
    # We need to reindex: table[L*4+C*2+R] -> f[L*1+C*2+R*4]
    reindexed = np.zeros(8, dtype=np.uint8)
    for L in range(2):
        for C in range(2):
            for R in range(2):
                old_idx = (L << 2) | (C << 1) | R  # ECA convention
                new_idx = L | (C << 1) | (R << 2)   # ANF convention
                reindexed[new_idx] = table[old_idx]
    
    anf = mobius_transform(reindexed)
    
    support = set()
    support_names = set()
    degree_counts = [0, 0, 0, 0]  # degrees 0, 1, 2, 3
    
    for idx in range(8):
        if anf[idx]:
            support.add(idx)
            support_names.add(MONOMIAL_NAMES[idx])
            degree_counts[MONOMIAL_DEGREES[idx]] += 1
    
    degree_profile = tuple(degree_counts)
    algebraic_degree = max((MONOMIAL_DEGREES[idx] for idx in support), default=-1)
    
    return {
        'rule': rule_number,
        'truth_table': table,
        'anf_coeffs': anf,
        'support': frozenset(support),
        'support_names': frozenset(support_names),
        'degree_profile': degree_profile,
        'algebraic_degree': algebraic_degree,
        'n_monomials': len(support),
    }


def degree_profile_str(dp: tuple) -> str:
    """Human-readable degree profile string."""
    return f"({dp[0]},{dp[1]},{dp[2]},{dp[3]})"


def anf_polynomial_str(anf_coeffs: np.ndarray) -> str:
    """Human-readable polynomial string for an ANF."""
    terms = []
    for idx in range(8):
        if anf_coeffs[idx]:
            terms.append(MONOMIAL_NAMES[idx])
    if not terms:
        return "0"
    return " ⊕ ".join(terms)


def compute_all_anfs() -> list:
    """Compute ANF for all 256 ECA rules."""
    return [compute_anf(r) for r in range(256)]


def verify_anf(rule_number: int) -> bool:
    """Verify that the ANF correctly reproduces the truth table.
    
    Evaluates the polynomial at all 8 input combinations and checks
    against the original rule lookup table.
    """
    info = compute_anf(rule_number)
    table = info['truth_table']
    anf = info['anf_coeffs']
    
    for L in range(2):
        for C in range(2):
            for R in range(2):
                # Evaluate polynomial at (L, C, R)
                result = 0
                for idx in range(8):
                    if anf[idx]:
                        # Check if all variables in this monomial are 1
                        # idx bits: bit0=L, bit1=C, bit2=R
                        vars_present = True
                        if idx & 0b001 and not L:
                            vars_present = False
                        if idx & 0b010 and not C:
                            vars_present = False
                        if idx & 0b100 and not R:
                            vars_present = False
                        if vars_present:
                            result ^= 1
                
                eca_idx = (L << 2) | (C << 1) | R
                if result != table[eca_idx]:
                    return False
    return True
