import random
from pathlib import Path

import pytest

from kersplit.berlekamp import factor
from kersplit.polyio import parse_problem
from tests.helpers import LARGE_PRIME, PRIME_512, as_pairs, random_factored

EXAMPLES = Path(__file__).parent.parent / "examples"
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 101]
CASES_PER_SMALL_PRIME_AND_METHOD = 25
LARGE_PRIME_CASES = {LARGE_PRIME: 130, PRIME_512: 20}


def _check(f, method, seed=0):
    ours = factor(f, method=method, seed=seed)
    sage = f.factor()
    assert ours.unit == sage.unit()
    assert as_pairs(ours.factors) == as_pairs(sage)
    assert ours.expand() == f


@pytest.mark.parametrize("p", SMALL_PRIMES)
@pytest.mark.parametrize("method", ["deterministic", "random"])
def test_random_polynomials_small_p(p, method):
    rng = random.Random(1000 * p + len(method))
    for _ in range(CASES_PER_SMALL_PRIME_AND_METHOD):
        f, _ = random_factored(p, rng, max_factors=6, max_degree=4, multiplicities=[1, 1, 2, 3, p, 2 * p])
        _check(f, method, seed=rng.randrange(10**6))


@pytest.mark.parametrize("p", sorted(LARGE_PRIME_CASES), ids=["987013", "512-bit"])
def test_random_polynomials_large_p(p):
    rng = random.Random(p % 10**6)
    for _ in range(LARGE_PRIME_CASES[p]):
        f, _ = random_factored(p, rng, max_factors=6, max_degree=4, multiplicities=[1, 1, 2, 3])
        _check(f, "random", seed=rng.randrange(10**6))


@pytest.mark.parametrize("sample", [1, 2, 3])
def test_hw3_samples(sample):
    problem = parse_problem((EXAMPLES / f"hw3_sample{sample}.txt").read_text())
    _check(problem.f, "random")


def test_multiplicities_p_2p_and_p_squared():
    rng = random.Random(7)
    f, _ = random_factored(13, rng, max_factors=3, max_degree=2, multiplicities=[13, 26, 169])
    _check(f, "deterministic")
    _check(f, "random")
