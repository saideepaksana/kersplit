from sage.rings.polynomial.polynomial_element import Polynomial


def pth_root(f: Polynomial) -> Polynomial:
    """Return h with h^p = f (notes, Claim 12): keep the coefficients of x^0, x^p, x^2p, ..."""
    p = int(f.base_ring().characteristic())
    coeffs = f.list()
    if any(c != 0 for i, c in enumerate(coeffs) if i % p):
        raise ValueError(f"{f} is not a p-th power")
    return f.parent()(coeffs[::p])


def squarefree_decomposition(f: Polynomial) -> list[tuple[Polynomial, int]]:
    """Split monic f into pairs (g, e) with g square-free, pairwise coprime and f = prod g^e.

    gcd(f, f') removes one copy of every factor whose multiplicity is not divisible by p.
    Factors whose multiplicity is divisible by p survive in c, which is then a p-th power.
    """
    if f.degree() < 1:
        return []
    p = int(f.base_ring().characteristic())
    c = f.gcd(f.derivative())
    w = f // c
    parts = []
    multiplicity = 1
    while w.degree() > 0:
        y = w.gcd(c)
        exact = w // y
        if exact.degree() > 0:
            parts.append((exact, multiplicity))
        multiplicity += 1
        w = y
        c = c // y
    if c.degree() > 0:
        parts += [(g, e * p) for g, e in squarefree_decomposition(pth_root(c))]
    return parts
