from pathlib import Path

import pytest

from kersplit.polyio import parse_problem

EXAMPLES = Path(__file__).parent.parent / "examples"


def test_parses_multiline_hw3_sample():
    problem = parse_problem((EXAMPLES / "hw3_sample2.txt").read_text())
    assert problem.p == 987013
    assert problem.f.degree() == 24
    assert problem.f.list()[0] == 805981


@pytest.mark.parametrize(
    "text, message",
    [
        ("p:10\nd:1\nf:x+1\n", "not a prime"),
        ("p:7\nd:2\nf:x+1\n", "degree"),
        ("p:7\nd:1\n", "missing"),
        ("p:7\nd:1\nf:x+y\n", "cannot read"),
        ("p:7\nd:0\nf:0\n", "non-zero"),
        ("p:7\nd:1\nf:x/7\n", "cannot read"),
        ("p:7\np:5\nd:1\nf:x\n", "more than once"),
    ],
)
def test_rejects_bad_input(text, message):
    with pytest.raises(ValueError, match=message):
        parse_problem(text)
