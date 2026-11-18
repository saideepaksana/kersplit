from sage.rings.polynomial.polynomial_element import Polynomial


def power_mod(g: Polynomial, e: int, f: Polynomial) -> Polynomial:
    """Return g^e mod f by repeated squaring."""
    result = f.parent().one()
    base = g % f
    e = int(e)
    while e > 0:
        if e & 1:
            result = result * base % f
        base = base * base % f
        e >>= 1
    return result
