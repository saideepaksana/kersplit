import random

import pytest
from sage.all import GF, matrix, vector

from kersplit.kernel import left_kernel_basis


@pytest.mark.parametrize("p", [2, 3, 101, 2**61 - 1])
def test_matches_sage_on_random_low_rank_matrices(p):
    rng = random.Random(p)
    field = GF(p)
    for _ in range(30):
        n = rng.randint(1, 8)
        rank = rng.randint(0, n)
        a = matrix(field, n, rank, [rng.randrange(p) for _ in range(n * rank)])
        b = matrix(field, rank, n, [rng.randrange(p) for _ in range(n * rank)])
        m = a * b if rank else matrix(field, n, n)
        basis = left_kernel_basis(m)
        assert len(basis) == m.left_kernel().dimension()
        for v in basis:
            assert vector(field, v) * m == 0
        if basis:
            assert matrix(field, basis).rank() == len(basis)


def test_identity_has_trivial_kernel():
    assert left_kernel_basis(matrix(GF(7), 3, 3, 1)) == []
