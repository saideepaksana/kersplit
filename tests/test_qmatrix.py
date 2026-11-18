import random

import pytest

from kersplit.qmatrix import berlekamp_matrix
from tests.helpers import ring


def test_x4_plus_1_mod_3():
    x = ring(3).gen()
    q = berlekamp_matrix(x**4 + 1)
    assert [list(row) for row in q] == [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 2, 0], [0, 1, 0, 0]]


@pytest.mark.parametrize("p", [2, 7, 987013])
def test_rows_are_powers_of_x(p):
    r = ring(p)
    x = r.gen()
    rng = random.Random(p)
    f = r([rng.randrange(p) for _ in range(9)] + [1])
    q = berlekamp_matrix(f)
    for i in range(f.degree()):
        assert r(list(q[i])) == pow(x, i * p, f)
