"""Commutator operators for ECA analysis.

DEFINITIONS (from Wet Math / Groovy Commutator):

  φ(S)       = Apply the CA rule table. The raw rules. Pure function.
  D(S)       = S ⊕ φ(S). "Which cells are about to change?" NOT spatial adjacency.
  I(S, d)    = S ⊕ d. "Apply the change pattern d to state S."
  E(S)       = I(S, D(S)) = S ⊕ (S ⊕ φ(S)) = φ(S). Evolution = integration of differentiation.
  C(a, b)    = a ⊕ b. Comparison (XOR for binary CAs).
  D₂(S)     = D(D(S)). Differentiate the derivative-as-state. NOT the second time derivative.
             = D(S) ⊕ φ(D(S)). Feed the change pattern into the rules as if it were a world.

  G(S) = C( D(E(S)), E(D(S)) )
       = C( D(φ(S)), I(D(S), D₂(S)) )
       = C( D(φ(S)), D(S) ⊕ D₂(S) )

  Path 1: D(E(S)) — evolve first, then ask what's changing in the evolved state.
  Path 2: E(D(S)) — differentiate first, then evolve the change pattern.

Note: In CAs, E(S) = φ(S) identically. But writing E = I(S, D(S)) matters because
Path 2 requires D₂ — the system dreaming about its own dynamics.
"""

import numpy as np
from src.eca import step


# ═══ Primitives ═══

def phi(state: np.ndarray, rule_number: int) -> np.ndarray:
    """φ(S): Apply the CA rule table. The raw rules."""
    return step(state, rule_number)


def D(state: np.ndarray, rule_number: int) -> np.ndarray:
    """D(S) = S ⊕ φ(S). Which cells are about to change?

    This is the Groovy Commutator's differentiation operator.
    NOT spatial adjacency — it detects which cells will flip on the next step.
    """
    return state ^ phi(state, rule_number)


def I(state: np.ndarray, delta: np.ndarray) -> np.ndarray:
    """I(S, d) = S ⊕ d. Apply change pattern d to state S."""
    return state ^ delta


def C(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """C(a, b) = a ⊕ b. Comparison operator."""
    return a ^ b


# ═══ Derived Operators ═══

def E(state: np.ndarray, rule_number: int) -> np.ndarray:
    """E(S) = I(S, D(S)) = φ(S). Evolution operator (composed from I and D)."""
    return phi(state, rule_number)


def E_n(state: np.ndarray, rule_number: int, n: int) -> np.ndarray:
    """E^n(S): Apply evolution n times."""
    s = state.copy()
    for _ in range(n):
        s = E(s, rule_number)
    return s


def D2(state: np.ndarray, rule_number: int) -> np.ndarray:
    """D₂(S) = D(D(S)). Differentiate the derivative-as-state.

    NOT the second time derivative!
    D₂ = D(S) ⊕ φ(D(S)). Feed the change pattern into φ as if it were a world.
    "The system dreaming about its own dynamics."
    """
    d_s = D(state, rule_number)
    return D(d_s, rule_number)  # = d_s ⊕ φ(d_s)


# ═══ The Groovy Commutator ═══

def G(state: np.ndarray, rule_number: int) -> np.ndarray:
    """G(S) = C( D(E(S)), E(D(S)) ).

    Path 1: D(E(S)) — evolve, then differentiate the result.
    Path 2: E(D(S)) = I(D(S), D₂(S)) — evolve the derivative-as-state.

    In CA terms:
      Path 1: D(φ(S)) = φ(S) ⊕ φ(φ(S))
      Path 2: I(D(S), D₂(S)) = D(S) ⊕ D₂(S)
    """
    evolved = E(state, rule_number)
    d_of_e = D(evolved, rule_number)          # Path 1: D(E(S))

    d_s = D(state, rule_number)
    e_of_d = E(d_s, rule_number)              # Path 2: E(D(S)) = φ(D(S))
    # Equivalently: I(d_s, D2(state, rule_number))

    return C(d_of_e, e_of_d)


# ═══ Commutator Variants ═══

def commutator_ED2(state: np.ndarray, rule_number: int) -> np.ndarray:
    """[E, D₂](S) = D₂(E(S)) ⊕ E(D₂(S))."""
    return C(D2(E(state, rule_number), rule_number),
             E(D2(state, rule_number), rule_number))


def commutator_EnD(state: np.ndarray, rule_number: int, n: int) -> np.ndarray:
    """[E^n, D](S) = D(E^n(S)) ⊕ E^n(D(S))."""
    en_s = E_n(state, rule_number, n)
    d_s = D(state, rule_number)
    en_d = E_n(d_s, rule_number, n)
    return C(D(en_s, rule_number), en_d)


def commutator_nested(state: np.ndarray, rule_number: int) -> np.ndarray:
    """[E, [E, D]](S) = [E,D](E(S)) ⊕ E([E,D](S)).

    The nested commutator: commutator of E with the commutator itself.
    """
    g_of_e = G(E(state, rule_number), rule_number)  # [E,D](E(S))
    e_of_g = E(G(state, rule_number), rule_number)  # E([E,D](S))
    return C(g_of_e, e_of_g)


def anti_commutator(state: np.ndarray, rule_number: int) -> np.ndarray:
    """{E,D}(S) = D(E(S)) AND E(D(S)).

    Where both paths agree on producing a 1.
    """
    evolved = E(state, rule_number)
    d_of_e = D(evolved, rule_number)
    d_s = D(state, rule_number)
    e_of_d = E(d_s, rule_number)
    return d_of_e & e_of_d


# ═══ Registry ═══
COMMUTATOR_VARIANTS = {
    "[E,D]": G,
    "[E,D₂]": commutator_ED2,
    "[E²,D]": lambda s, r: commutator_EnD(s, r, 2),
    "[E³,D]": lambda s, r: commutator_EnD(s, r, 3),
    "[E⁴,D]": lambda s, r: commutator_EnD(s, r, 4),
    "[E,[E,D]]": commutator_nested,
    "{E,D}": anti_commutator,
}
