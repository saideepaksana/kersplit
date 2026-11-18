import re
from dataclasses import dataclass

from sage.all import GF, PolynomialRing, is_prime
from sage.rings.polynomial.polynomial_element import Polynomial

_FIELD = re.compile(r"^\s*([pdf])\s*:", re.MULTILINE)


@dataclass(frozen=True)
class Problem:
    p: int
    f: Polynomial


def parse_problem(text: str) -> Problem:
    """Parse HW3-format input: lines 'p:', 'd:' and 'f:', where f may continue over several lines."""
    fields = _split_fields(text)
    missing = {"p", "d", "f"} - fields.keys()
    if missing:
        raise ValueError(f"missing field(s): {', '.join(sorted(missing))}")
    p = _parse_int(fields["p"], "p")
    d = _parse_int(fields["d"], "d")
    if p < 2 or not is_prime(p):
        raise ValueError(f"p = {p} is not a prime")
    f = parse_polynomial(fields["f"], p)
    if f.is_zero():
        raise ValueError("f must be non-zero")
    if f.degree() != d:
        raise ValueError(f"f has degree {f.degree()}, but d = {d}")
    return Problem(p, f)


def parse_polynomial(text: str, p: int) -> Polynomial:
    ring = PolynomialRing(GF(p), "x")
    try:
        return ring(" ".join(text.split()))
    except (TypeError, SyntaxError, NameError, ArithmeticError) as error:
        raise ValueError(f"cannot read polynomial {text!r}") from error


def _split_fields(text: str) -> dict[str, str]:
    matches = list(_FIELD.finditer(text))
    names = [m.group(1) for m in matches]
    repeated = sorted({name for name in names if names.count(name) > 1})
    if repeated:
        raise ValueError(f"field(s) given more than once: {', '.join(repeated)}")
    ends = [m.start() for m in matches[1:]] + [len(text)]
    return {m.group(1): text[m.end() : end].strip() for m, end in zip(matches, ends)}


def _parse_int(text: str, name: str) -> int:
    try:
        return int(text)
    except ValueError as error:
        raise ValueError(f"{name} must be an integer, got {text!r}") from error
