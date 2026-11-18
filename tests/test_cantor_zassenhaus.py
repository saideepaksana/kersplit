import random

import pytest

from bench.cantor_zassenhaus import cantor_zassenhaus
from tests.helpers import LARGE_PRIME, as_pairs, random_factored


@pytest.mark.parametrize("p", [2, 3, 101, LARGE_PRIME])
def test_baseline_matches_sage(p):
    rng = random.Random(p)
    for _ in range(20):
        f, known = random_factored(p, rng, max_factors=6, max_degree=4, multiplicities=[1])
        assert as_pairs((h, 1) for h in cantor_zassenhaus(f.monic(), rng)) == as_pairs(known)
