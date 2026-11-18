from sage.all import matrix
from sage.matrix.matrix0 import Matrix
from sage.rings.polynomial.polynomial_element import Polynomial

from kersplit.polyops import power_mod


def berlekamp_matrix(f: Polynomial) -> Matrix:
    """Return the d x d matrix Q whose row i holds the coefficients of x^(ip) mod f."""
    d = f.degree()
    p = int(f.base_ring().characteristic())
    x_p = power_mod(f.parent().gen(), p, f)
    rows = []
    row = f.parent().one()
    for _ in range(d):
        rows.append(row.padded_list(d))
        row = row * x_p % f
    return matrix(f.base_ring(), rows)
