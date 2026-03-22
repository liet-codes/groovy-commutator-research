"""Affine decomposition of ECA rules (Brooklyn Rose approach).

Each ECA rule is a Boolean function f: {0,1}^3 -> {0,1}.
We decompose f = a + p (mod 2) where:
  a = best-fit affine (XOR) function
  p = residual perturbation

The 8 affine functions on 3 variables are all XOR combinations:
  0, L, C, R, L^C, L^R, C^R, L^C^R  (plus complements = 16 total affine functions)
"""

import numpy as np
from src.eca import rule_lookup


# All 16 affine Boolean functions on 3 inputs (including complements)
# Inputs indexed as: [L*4 + C*2 + R] for the 8 neighborhoods
def _build_affine_basis():
    """Build all 16 affine functions on 3 binary inputs."""
    neighborhoods = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
                               [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]], dtype=np.uint8)
    L, C, R = neighborhoods[:, 0], neighborhoods[:, 1], neighborhoods[:, 2]

    affine_funcs = {}
    # All XOR subsets of {L, C, R} with optional complement
    subsets = [
        ("0",       np.zeros(8, dtype=np.uint8)),
        ("L",       L),
        ("C",       C),
        ("R",       R),
        ("L^C",     L ^ C),
        ("L^R",     L ^ R),
        ("C^R",     C ^ R),
        ("L^C^R",   L ^ C ^ R),
    ]
    for name, func in subsets:
        affine_funcs[name] = func
        affine_funcs["1^" + name if name != "0" else "1"] = 1 - func

    return affine_funcs


AFFINE_FUNCTIONS = _build_affine_basis()


def affine_decompose(rule_number: int):
    """Decompose an ECA rule into best-fit affine + perturbation.

    Returns dict with:
      - 'rule': rule number
      - 'table': the 8-entry rule lookup table
      - 'best_affine_name': name of the closest affine function
      - 'best_affine_table': the affine function's table
      - 'perturbation': XOR difference (residual)
      - 'perturbation_magnitude': number of differing entries (0-8)
      - 'perturbation_weight': perturbation_magnitude / 8
      - 'is_affine': whether the rule is exactly affine
    """
    table = rule_lookup(rule_number)

    best_name = None
    best_dist = 9
    best_affine = None

    for name, afunc in AFFINE_FUNCTIONS.items():
        dist = int(np.sum(table != afunc))
        if dist < best_dist:
            best_dist = dist
            best_name = name
            best_affine = afunc

    perturbation = table ^ best_affine

    return {
        'rule': rule_number,
        'table': table,
        'best_affine_name': best_name,
        'best_affine_table': best_affine.copy(),
        'perturbation': perturbation,
        'perturbation_magnitude': best_dist,
        'perturbation_weight': best_dist / 8.0,
        'is_affine': best_dist == 0,
    }


def perturbation_coherence(perturbation: np.ndarray) -> float:
    """Measure how "structured" the perturbation is.

    Returns a coherence score based on how localized the perturbation is
    in the input space. A perturbation affecting adjacent neighborhoods
    scores higher than one affecting scattered neighborhoods.

    Score: 1 - (number of transitions in perturbation pattern / max transitions).
    """
    if np.sum(perturbation) == 0:
        return 0.0
    transitions = int(np.sum(perturbation[:-1] != perturbation[1:]))
    max_transitions = len(perturbation) - 1
    if max_transitions == 0:
        return 1.0
    return 1.0 - transitions / max_transitions


def classify_all_rules():
    """Decompose all 256 ECA rules and return a list of results."""
    results = []
    for r in range(256):
        dec = affine_decompose(r)
        dec['perturbation_coherence'] = perturbation_coherence(dec['perturbation'])
        results.append(dec)
    return results
