import random
from collections.abc import Callable

from sage.rings.polynomial.polynomial_element import Polynomial

from kersplit.polyops import power_mod


def split_random(
    f: Polynomial,
    basis: list[Polynomial],
    rng: random.Random,
    log: Callable[[str], None],
) -> list[Polynomial]:
    """Return the monic irreducible factors of square-free f, using random elements of B.

    A random g in B is a random constant c_i modulo each irreducible h_i. For odd p,
    g^((p-1)/2) is 1 mod h_i exactly when c_i is a non-zero square, which happens with
    probability (p-1)/(2p), independently for each i. So gcd(h, g^((p-1)/2) - 1) splits any
    reducible h with probability close to 1/2 per round. For p = 2, c_i is 0 or 1, and gcd(h, g)
    splits the same way.
    """
    p = int(f.base_ring().characteristic())
    field = f.base_ring()
    factors = [f]
    rounds = 0
    while len(factors) < len(basis):
        rounds += 1
        g = sum(field(rng.randrange(p)) * b for b in basis)
        factors = [part for h in factors for part in split_once(h, g, p)]
        log(f"  round {rounds}: {len(factors)} of {len(basis)} factors found")
    return factors


def split_once(h: Polynomial, g: Polynomial, p: int) -> list[Polynomial]:
    """Return [d, h // d] for d = gcd(h, g^((p-1)/2) - 1) (gcd(h, g) if p = 2) if d is proper, else [h]."""
    if h.degree() <= 1:
        return [h]
    a = g % h
    candidate = a if p == 2 else power_mod(a, (p - 1) // 2, h) - 1
    d = h.gcd(candidate)
    if 0 < d.degree() < h.degree():
        return [d, h // d]
    return [h]
