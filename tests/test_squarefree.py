import random

import pytest

from kersplit.squarefree import pth_root, squarefree_decomposition
from tests.helpers import random_factored, ring


def test_pth_root_inverts_frobenius():
    r = ring(5)
    x = r.gen()
    h = x**3 + 2 * x + 4
    assert pth_root(h**5) == h


def test_pth_root_rejects_non_power():
    x = ring(5).gen()
    with pytest.raises(ValueError):
        pth_root(x**5 + x)


def test_derivative_zero_case():
    x = ring(3).gen()
    f = (x + 1) ** 3 * (x**2 + 1) ** 6
    assert squarefree_decomposition(f) == [(x + 1, 3), (x**2 + 1, 6)]


@pytest.mark.parametrize("p", [2, 3, 5, 7, 101])
def test_matches_sage_squarefree_decomposition(p):
    rng = random.Random(p)
    for _ in range(40):
        f, _ = random_factored(p, rng, max_factors=4, max_degree=3, multiplicities=[1, 2, 3, p, p + 1, 2 * p])
        ours = sorted((tuple(g.list()), e) for g, e in squarefree_decomposition(f.monic()))
        sage = sorted((tuple(g.list()), int(e)) for g, e in f.monic().squarefree_decomposition())
        assert ours == sage
