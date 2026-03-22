"""Quantitative measures for binary state arrays."""

import gzip
import numpy as np
from scipy import signal


def shannon_entropy(state: np.ndarray) -> float:
    """Shannon entropy of a binary array (bits).

    Returns entropy in [0, 1] — 0 for uniform, 1 for 50/50.
    """
    n = len(state)
    if n == 0:
        return 0.0
    p1 = np.sum(state) / n
    p0 = 1.0 - p1
    if p0 <= 0 or p1 <= 0:
        return 0.0
    return -(p0 * np.log2(p0) + p1 * np.log2(p1))


def block_entropy(state: np.ndarray, block_size: int = 4) -> float:
    """Shannon entropy of non-overlapping blocks of given size.

    Captures spatial structure beyond single-cell density.
    """
    n = len(state)
    nblocks = n // block_size
    if nblocks == 0:
        return 0.0
    trimmed = state[:nblocks * block_size].reshape(nblocks, block_size)
    # Convert each block to an integer
    powers = 1 << np.arange(block_size)
    block_ids = trimmed @ powers
    # Count frequencies
    unique, counts = np.unique(block_ids, return_counts=True)
    probs = counts / counts.sum()
    return -np.sum(probs * np.log2(probs))


def spatial_correlation_length(state: np.ndarray) -> float:
    """Spatial correlation length via autocorrelation decay.

    Returns the 1/e decay distance of the spatial autocorrelation function.
    """
    n = len(state)
    if n < 4:
        return 0.0
    centered = state.astype(np.float64) - np.mean(state)
    var = np.var(state)
    if var < 1e-12:
        return 0.0
    # Full autocorrelation via FFT
    acf = np.real(np.fft.ifft(np.abs(np.fft.fft(centered)) ** 2)) / (n * var)
    # Find first crossing below 1/e
    half = n // 2
    acf_half = acf[1:half]  # skip lag 0
    threshold = 1.0 / np.e
    below = np.where(acf_half < threshold)[0]
    if len(below) == 0:
        return float(half)
    return float(below[0] + 1)


def temporal_autocorrelation(time_series: np.ndarray, max_lag: int = 50) -> np.ndarray:
    """Temporal autocorrelation of a sequence of scalar values.

    Returns autocorrelation for lags 0..max_lag.
    """
    n = len(time_series)
    if n < 2:
        return np.array([1.0])
    centered = time_series - np.mean(time_series)
    var = np.var(time_series)
    if var < 1e-12:
        return np.ones(min(max_lag + 1, n))
    max_lag = min(max_lag, n - 1)
    acf = np.correlate(centered, centered, mode='full')
    acf = acf[n - 1:]  # keep non-negative lags
    acf = acf / acf[0]
    return acf[:max_lag + 1]


def temporal_autocorrelation_mean(time_series: np.ndarray, max_lag: int = 50) -> float:
    """Mean temporal autocorrelation (averaged over lags 1..max_lag)."""
    acf = temporal_autocorrelation(time_series, max_lag)
    if len(acf) < 2:
        return 0.0
    return float(np.mean(np.abs(acf[1:])))


def hamming_distance(a: np.ndarray, b: np.ndarray) -> int:
    """Hamming distance between two binary arrays."""
    return int(np.sum(a != b))


def hamming_distance_normalized(a: np.ndarray, b: np.ndarray) -> float:
    """Normalized Hamming distance in [0, 1]."""
    return hamming_distance(a, b) / len(a)


def compressibility_ratio(state: np.ndarray) -> float:
    """Compressibility ratio as proxy for Kolmogorov complexity.

    Returns compressed_size / original_size. Lower = more compressible = simpler.
    """
    raw = state.tobytes()
    if len(raw) == 0:
        return 1.0
    compressed = gzip.compress(raw, compresslevel=9)
    return len(compressed) / len(raw)
