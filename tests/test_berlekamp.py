import pytest

from kersplit.berlekamp import berlekamp_basis, factor
from tests.helpers import ring


@pytest.mark.parametrize("p", [3, 5, 7, 11, 13])
def test_x4_plus_1_has_two_factors_mod_every_odd_prime(p):
    x = ring(p).gen()
    assert len(berlekamp_basis(x**4 + 1)) == 2


def test_worked_example_basis():
    x = ring(3).gen()
    assert berlekamp_basis(x**4 + 1) == [1, x**3 + x]


def test_irreducible_has_one_dimensional_kernel():
    x = ring(2).gen()
    assert berlekamp_basis(x**5 + x**2 + 1) == [1]


def test_constant_and_linear():
    r = ring(7)
    x = r.gen()
    assert factor(r(3)).factors == []
    assert factor(3 * x + 1).factors == [(x + 5, 1)]


def test_unit_and_multiplicities():
    x = ring(5).gen()
    f = 2 * (x + 1) ** 5 * (x**2 + 2) ** 2 * (x + 3)
    result = factor(f)
    assert result.unit == 2
    assert result.factors == [(x + 1, 5), (x + 3, 1), (x**2 + 2, 2)]
    assert result.expand() == f


def test_rejects_non_prime_field():
    from sage.all import GF, PolynomialRing

    x = PolynomialRing(GF(4, "a"), "x").gen()
    with pytest.raises(ValueError, match="prime"):
        factor(x**2 + x + 1)


def test_rejects_zero_and_unknown_method():
    r = ring(5)
    with pytest.raises(ValueError):
        factor(r(0))
    with pytest.raises(ValueError):
        factor(r.gen(), method="magic")


@pytest.mark.parametrize("make", ["rational", "multivariate", "constant_element"])
def test_rejects_inputs_outside_univariate_z_p(make):
    from sage.all import GF, QQ, PolynomialRing

    polys = {
        "rational": lambda: PolynomialRing(QQ, "x").gen() ** 2 - 1,
        "multivariate": lambda: PolynomialRing(GF(5), "x,y").gen(0) ** 2 - 1,
        "constant_element": lambda: GF(5)(3),
    }
    with pytest.raises(ValueError):
        factor(polys[make]())


def test_deterministic_refused_for_huge_p():
    x = ring(1000003).gen()
    with pytest.raises(ValueError, match="random"):
        factor((x + 1) * (x + 2), method="deterministic")


def test_constant_expands_to_polynomial():
    r = ring(7)
    assert factor(r(3)).expand().parent() is r


def test_str_matches_sage_format():
    x = ring(2).gen()
    f = x**3 * (x**2 + x + 1) ** 5 * (x**3 + x + 1)
    assert str(factor(f)) == str(f.factor())
