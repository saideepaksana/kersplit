import random
from collections.abc import Callable
from dataclasses import dataclass

from sage.all import identity_matrix
from sage.matrix.matrix0 import Matrix
from sage.rings.polynomial.polynomial_element import Polynomial

from kersplit.kernel import left_kernel_basis
from kersplit.qmatrix import berlekamp_matrix
from kersplit.random_split import split_random
from kersplit.squarefree import squarefree_decomposition

DETERMINISTIC_MAX_P = 1000
DETERMINISTIC_LIMIT_P = 10**6
MAX_PRINTED_Q_DEGREE = 8
METHODS = ("auto", "deterministic", "random")


@dataclass(frozen=True)
class Factorization:
    """f = unit * prod(h^e for h, e in factors), with distinct monic irreducible h."""

    unit: object
    factors: list[tuple[Polynomial, int]]

    def expand(self) -> Polynomial:
        product = self.unit
        for h, e in self.factors:
            product *= h**e
        return product

    def __str__(self) -> str:
        terms = [_power(h, e) for h, e in self.factors]
        if self.unit != 1 or not terms:
            terms.insert(0, str(self.unit))
        return " * ".join(terms)


def berlekamp_basis(f: Polynomial, q: Matrix | None = None) -> list[Polynomial]:
    """Return a basis of B = {g : g^p = g mod f}, the kernel of Q - I; basis[0] = 1."""
    if q is None:
        q = berlekamp_matrix(f)
    vectors = left_kernel_basis(q - identity_matrix(f.base_ring(), f.degree()))
    return [f.parent()(v) for v in vectors]


def split_deterministic(
    f: Polynomial, basis: list[Polynomial], log: Callable[[str], None]
) -> list[Polynomial]:
    """Return the monic irreducible factors of square-free f.

    g^p - g = prod_{s in Z_p} (g - s), so for g in B every factor h of f is the product of the
    pairwise coprime gcd(h, g - s). Some basis element separates any two irreducible factors,
    so refining by every basis element yields all r of them.
    """
    p = int(f.base_ring().characteristic())
    factors = [f]
    for g in basis[1:]:
        if len(factors) == len(basis):
            break
        log(f"  split with g = {g}")
        factors = [part for h in factors for part in _split_by(h, g, p, log)]
    return factors


def _split_by(h: Polynomial, g: Polynomial, p: int, log: Callable[[str], None]) -> list[Polynomial]:
    if h.degree() <= 1:
        return [h]
    g = g % h
    log(f"    piece {h}:")
    parts = []
    remaining = h.degree()
    for s in range(p):
        if remaining == 0:
            break
        d = h.gcd(g - s)
        if d.degree() > 0:
            log(f"      s = {s}: gcd = {d}")
            parts.append(d)
            remaining -= d.degree()
    return parts


def factor(
    f: Polynomial, method: str = "auto", trace: bool = False, seed: int = 0
) -> Factorization:
    """Factor a non-zero f in Z_p[x] completely with Berlekamp's algorithm."""
    if not isinstance(f, Polynomial):
        raise ValueError("f must be a univariate polynomial")
    field = f.base_ring()
    if not (field.is_prime_field() and field.characteristic() > 0):
        raise ValueError("f must have coefficients in Z_p for a prime p")
    if f.is_zero():
        raise ValueError("cannot factor the zero polynomial")
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}")
    p = int(field.characteristic())
    if method == "deterministic" and p > DETERMINISTIC_LIMIT_P:
        raise ValueError(
            f"the deterministic split tries all p values of s; use method='random' for p > {DETERMINISTIC_LIMIT_P}"
        )
    log = print if trace else _quiet
    use_random = method == "random" or (method == "auto" and p > DETERMINISTIC_MAX_P)
    rng = random.Random(seed)
    log(f"f = {f} over Z_{p}, degree {f.degree()}")

    factors = []
    for part, e in squarefree_decomposition(f.monic()):
        log(f"square-free part with multiplicity {e}: {part}")
        factors += [(h, e) for h in _factor_squarefree(part, use_random, rng, log)]
    factors.sort(key=lambda he: (he[0].degree(), he[0].list()))
    result = Factorization(f.parent()(f.leading_coefficient()), factors)
    log(f"result: {result}")
    return result


def _factor_squarefree(
    f: Polynomial, use_random: bool, rng: random.Random, log: Callable[[str], None]
) -> list[Polynomial]:
    if f.degree() == 1:
        return [f]
    q = berlekamp_matrix(f)
    if f.degree() <= MAX_PRINTED_Q_DEGREE:
        log(f"  Q (row i = x^(ip) mod f):\n{_indent(str(q))}")
    else:
        log(f"  Q is {f.degree()} x {f.degree()} (row i = x^(ip) mod f)")
    basis = berlekamp_basis(f, q)
    log(f"  kernel of Q - I: basis {{{', '.join(map(str, basis))}}} -> r = {len(basis)}")
    if len(basis) == 1:
        log("  r = 1, so this part is irreducible")
        return [f]
    if use_random:
        return split_random(f, basis, rng, log)
    return split_deterministic(f, basis, log)


def _power(h: Polynomial, e: int) -> str:
    base = str(h) if len(h.exponents()) == 1 else f"({h})"
    return base if e == 1 else f"{base}^{e}"


def _indent(text: str) -> str:
    return "\n".join("    " + line for line in text.splitlines())


def _quiet(_: str) -> None:
    pass
