"""Commutator operators for ECA analysis."""

import numpy as np
from src.eca import step


def D(state: np.ndarray) -> np.ndarray:
    """Differentiation: XOR of adjacent cells (spatial derivative).

    D(S)[i] = S[i] XOR S[i+1], with periodic boundary.
    """
    return state ^ np.roll(state, -1)


def D2(state: np.ndarray) -> np.ndarray:
    """Second-order differentiation: D applied twice."""
    return D(D(state))


def I(state: np.ndarray) -> np.ndarray:
    """Integration: cumulative XOR (spatial integral).

    I(S)[i] = S[0] XOR S[1] XOR ... XOR S[i].
    """
    result = np.zeros_like(state)
    acc = np.uint8(0)
    for i in range(len(state)):
        acc ^= state[i]
        result[i] = acc
    return result


def E(state: np.ndarray, rule_number: int) -> np.ndarray:
    """Evolution: one CA step."""
    return step(state, rule_number)


def E_n(state: np.ndarray, rule_number: int, n: int) -> np.ndarray:
    """Multi-step evolution: apply E n times."""
    s = state.copy()
    for _ in range(n):
        s = step(s, rule_number)
    return s


def G(state: np.ndarray, rule_number: int) -> np.ndarray:
    """Groovy Commutator: G(S) = D(E(S)) XOR E(D(S))."""
    de = D(E(state, rule_number))
    ed = E(D(state), rule_number)
    return de ^ ed


def commutator_ED2(state: np.ndarray, rule_number: int) -> np.ndarray:
    """Second-order commutator: [E, D²](S) = D²(E(S)) XOR E(D²(S))."""
    return D2(E(state, rule_number)) ^ E(D2(state), rule_number)


def commutator_EI(state: np.ndarray, rule_number: int) -> np.ndarray:
    """Integration commutator: [E, I](S) = I(E(S)) XOR E(I(S))."""
    return I(E(state, rule_number)) ^ E(I(state), rule_number)


def commutator_EnD(state: np.ndarray, rule_number: int, n: int) -> np.ndarray:
    """Multi-step commutator: [E^n, D](S) = D(E^n(S)) XOR E^n(D(S))."""
    return D(E_n(state, rule_number, n)) ^ E_n(D(state), rule_number, n)


def commutator_nested(state: np.ndarray, rule_number: int) -> np.ndarray:
    """Nested commutator: [E, [E, D]](S).

    Inner: [E,D](S) = G(S)
    Outer: [E, G](S) = G(E(S)) XOR E(G(S))
    Wait -- need to be careful. [E, [E,D]] means we treat [E,D] as an operator.
    [E, F](S) = F(E(S)) XOR E(F(S)) where F = [E,D] = G
    """
    inner = G(state, rule_number)  # [E,D](S)
    # [E, [E,D]](S) = [E,D](E(S)) XOR E([E,D](S))
    g_of_e = G(E(state, rule_number), rule_number)
    e_of_g = E(inner, rule_number)
    return g_of_e ^ e_of_g


def anti_commutator(state: np.ndarray, rule_number: int) -> np.ndarray:
    """Anti-commutator {E,D}: Hamming-sum-based.

    {E,D}(S) = D(E(S)) AND E(D(S))  (bitwise AND as "sum" analogue).
    We use AND since in GF(2), the anti-commutator analog is the AND of both branches,
    capturing where both orderings agree on producing a 1.
    """
    de = D(E(state, rule_number))
    ed = E(D(state), rule_number)
    return de & ed


# Registry for easy iteration
COMMUTATOR_VARIANTS = {
    "[E,D]": G,
    "[E,D²]": commutator_ED2,
    "[E,I]": commutator_EI,
    "[E²,D]": lambda s, r: commutator_EnD(s, r, 2),
    "[E³,D]": lambda s, r: commutator_EnD(s, r, 3),
    "[E⁴,D]": lambda s, r: commutator_EnD(s, r, 4),
    "[E,[E,D]]": commutator_nested,
    "{E,D}": anti_commutator,
}
