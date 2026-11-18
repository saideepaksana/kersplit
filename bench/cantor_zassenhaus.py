"""Cantor-Zassenhaus for monic square-free f (notes, Ch. 7), ported from our Programming HW3 as a baseline.

Added over HW3: distinct-degree factorization stops once 2k > deg f, since the rest is irreducible.
"""

import random

from sage.rings.polynomial.polynomial_element import Polynomial

from kersplit.polyops import power_mod


def cantor_zassenhaus(f: Polynomial, rng: random.Random) -> list[Polynomial]:
    p = int(f.base_ring().characteristic())
    return [h for part, k in _distinct_degree(f, p) for h in _equal_degree(part, k, p, rng)]


def _distinct_degree(f: Polynomial, p: int) -> list[tuple[Polynomial, int]]:
    x = f.parent().gen()
    parts = []
    h = x
    k = 1
    while f.degree() > 0:
        if 2 * k > f.degree():
            parts.append((f, f.degree()))
            break
        h = power_mod(h, p, f)
        part = f.gcd(h - x)
        if part.degree() > 0:
            parts.append((part, k))
            f = f // part
            h = h % f
        k += 1
    return parts


def _equal_degree(f: Polynomial, k: int, p: int, rng: random.Random) -> list[Polynomial]:
    if f.degree() == k:
        return [f]
    ring = f.parent()
    while True:
        g = ring([rng.randrange(p) for _ in range(f.degree())])
        if p == 2:
            t = g % f
            s = t
            for _ in range(k - 1):
                t = t * t % f
                s += t
            d = f.gcd(s)
        else:
            d = f.gcd(power_mod(g, (p**k - 1) // 2, f) - 1)
        if 0 < d.degree() < f.degree():
            return _equal_degree(d, k, p, rng) + _equal_degree(f // d, k, p, rng)
