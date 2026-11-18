import random

from sage.all import GF, PolynomialRing
from sage.rings.polynomial.polynomial_element import Polynomial

LARGE_PRIME = 987013
PRIME_512 = 2**512 - 569


def ring(p: int):
    return PolynomialRing(GF(p), "x")


def random_monic_irreducible(p: int, degree: int, rng: random.Random) -> Polynomial:
    r = ring(p)
    while True:
        coeffs = [rng.randrange(p) for _ in range(degree)] + [1]
        h = r(coeffs)
        if h.is_irreducible():
            return h


def random_factored(
    p: int, rng: random.Random, max_factors: int, max_degree: int, multiplicities: list[int]
) -> tuple[Polynomial, list[tuple[Polynomial, int]]]:
    """Return f and its known factorization into distinct monic irreducibles."""
    known: dict[tuple, tuple[Polynomial, int]] = {}
    for _ in range(rng.randint(1, max_factors)):
        h = random_monic_irreducible(p, rng.randint(1, max_degree), rng)
        known.setdefault(tuple(h.list()), (h, rng.choice(multiplicities)))
    unit = GF(p)(rng.randrange(1, p))
    f = ring(p)(unit)
    for h, e in known.values():
        f *= h**e
    return f, sorted(known.values(), key=lambda he: (he[0].degree(), he[0].list()))


def as_pairs(factorization) -> list[tuple[tuple, int]]:
    return sorted((tuple(h.list()), int(e)) for h, e in factorization)
