# KerSplit design note

Written 2026-11-18. The source of truth for module boundaries and function signatures. Change it
in the same edit as the code it describes.

## Claim

`kersplit.factor(f)` returns the complete factorization of any non-zero `f` in `Z_p[x]`
(leading coefficient, distinct monic irreducibles, multiplicities) using Berlekamp's algorithm.
Sage supplies only the ring `GF(p)[x]` and its `+ - * % // gcd`; every algorithmic step is ours.

## Layout and ownership

```
kersplit/                 package, one module per step
  polyops.py      M1      power_mod (repeated squaring mod f)
  squarefree.py   M1      pth_root, squarefree_decomposition
  qmatrix.py      M1      berlekamp_matrix
  polyio.py       M1      parse_problem, parse_polynomial
  kernel.py       M2      left_kernel_basis (Gaussian elimination over Z_p)
  berlekamp.py    M2      berlekamp_basis, split_deterministic, factor
  random_split.py M3      split_random
  __main__.py     M2      command line: python -m kersplit
bench/            M3      cantor_zassenhaus.py (HW3 baseline), bench.py, results/
tests/            all     each member tests their own modules; M3 owns the cross-check
examples/                 HW3-format input files
demo.ipynb                presentation notebook; imports kersplit, contains no algorithm code
```

Dependencies flow one way: `polyops` <- `qmatrix`, `random_split` <- `berlekamp` (which also
uses `squarefree` and `kernel`) <- `__main__`, `bench`, `demo.ipynb`. `squarefree` and `kernel`
depend on nothing in the package.

## Interfaces

`Poly` is a Sage polynomial over `GF(p)`. A *square-free part* is monic with no repeated factor.

| Function | Contract |
|---|---|
| `power_mod(g, e, f) -> Poly` | `g^e mod f`, `e >= 0` |
| `pth_root(f) -> Poly` | `h` with `h^p = f`; raises `ValueError` if `f` is not a `p`-th power (notes, Claim 12) |
| `squarefree_decomposition(f) -> list[tuple[Poly, int]]` | for monic `f`: pairs `(g, e)`, `g` square-free, pairwise coprime, distinct `e`, `f = prod g^e` |
| `berlekamp_matrix(f) -> Matrix` | `d x d`, row `i` = coefficients of `x^(ip) mod f` |
| `left_kernel_basis(M) -> list[list]` | basis of `{v : v M = 0}`; first vector is `e_0` when row 0 of `M` is zero |
| `berlekamp_basis(f) -> list[Poly]` | basis of `{g : g^p = g mod f}`, `basis[0] = 1`; `len = r` |
| `split_deterministic(f, basis, log) -> list[Poly]` | all `r` monic irreducible factors of square-free `f`, via `gcd(h, g - s)` for `s` in `Z_p` |
| `split_once(h, g, p) -> list[Poly]` | `[d, h/d]` for `d = gcd(h, g^((p-1)/2) - 1)` (`gcd(h, g)` when `p = 2`) if `d` is a proper factor, else `[h]` |
| `split_random(f, basis, rng, log) -> list[Poly]` | same output as `split_deterministic`, by `split_once` with random `g` in `B` until there are `r` factors |
| `factor(f, method, trace, seed) -> Factorization` | `method` in `{"auto", "deterministic", "random"}`; `auto` is deterministic for `p <= 1000`; `deterministic` is refused for `p > 10^6`; raises `ValueError` unless `f` is a non-zero univariate polynomial over a prime field `Z_p` |

`log` is a `Callable[[str], None]`; `factor(trace=True)` passes `print`, otherwise a no-op.

## Invariants (each has a test)

1. `unit * prod h^e == f` and every `h` is irreducible (checked against Sage `factor()`).
2. `len(berlekamp_basis(f)) ==` number of distinct irreducible factors of square-free `f`.
3. Every kernel vector `v` satisfies `v (Q - I) = 0`, and the basis is linearly independent.
4. Output is deterministic: factors sorted by `(degree, coefficients)`; `split_random` uses a seeded
   `random.Random`.

## Acceptance

1. `pytest` passes: unit tests, the three HW3 samples, 500 random cross-checks against `factor()`
   for `p` in {2, 3, 5, 7, 11, 13, 101, 987013} and a 512-bit prime; multiplicities include `p` and
   `2p` for the small primes.
2. `python -m kersplit --trace < examples/x4plus1_mod3.txt` prints `Q`, the kernel basis
   `{1, x^3 + x}`, and `(x^2 + x + 2) * (x^2 + 2*x + 2)`.
3. `python -m bench.bench` writes CSV tables and plots to `bench/results/` in a few minutes.

## Alternatives considered

- **Sage `matrix.left_kernel()`**: one line, but it would hide the step the project is about.
  Kept only as a test oracle.
- **Our own `F_p[x]` arithmetic**: about 120 lines unrelated to Berlekamp, slower than FLINT, and
  it only duplicates HW1 to HW3 work.

## Out of scope

Factoring over Z or `F_(p^k)`, fast variants (Kaltofen-Shoup, Wiedemann), applications to codes.
