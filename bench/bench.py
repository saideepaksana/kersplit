"""Benchmarks for the talk. Writes CSV tables, summary.md and plots to bench/results/.

Run from the repository root:  python -m bench.bench
"""

import csv
import random
import statistics
import time
from collections.abc import Callable
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sage.all import GF, PolynomialRing, previous_prime
from sage.rings.polynomial.polynomial_element import Polynomial

from bench.cantor_zassenhaus import cantor_zassenhaus
from kersplit.berlekamp import berlekamp_basis, split_deterministic
from kersplit.qmatrix import berlekamp_matrix
from kersplit.random_split import split_random

RESULTS = Path(__file__).parent / "results"
TRIALS = 7
LARGE_PRIME = 987013
DEGREES = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
FIXED_DEGREE = 30
DETERMINISTIC_BITS = [2, 4, 6, 8, 10, 12, 14, 16]
RANDOMIZED_BITS = DETERMINISTIC_BITS + [32, 64, 128, 256, 512]

COLORS = {"berlekamp_random": "#2a78d6", "berlekamp_deterministic": "#eb6834",
          "cantor_zassenhaus": "#1baf7a", "sage_factor": "#eda100"}
LABELS = {"berlekamp_random": "Berlekamp (randomized)", "berlekamp_deterministic": "Berlekamp (deterministic)",
          "cantor_zassenhaus": "Cantor-Zassenhaus (HW3)", "sage_factor": "Sage factor() (C library)"}
STAGE_COLORS = {"q_matrix": "#2a78d6", "kernel": "#eb6834", "split": "#1baf7a"}
INK, MUTED, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"

Splitter = Callable[[Polynomial, random.Random], list[Polynomial]]


def berlekamp(f: Polynomial, rng: random.Random, deterministic: bool) -> tuple[list[Polynomial], dict]:
    """Factor square-free monic f, timing the three stages separately."""
    t0 = time.perf_counter()
    q = berlekamp_matrix(f)
    t1 = time.perf_counter()
    basis = berlekamp_basis(f, q)
    t2 = time.perf_counter()
    if len(basis) == 1:
        factors = [f]
    elif deterministic:
        factors = split_deterministic(f, basis, _quiet)
    else:
        factors = split_random(f, basis, rng, _quiet)
    t3 = time.perf_counter()
    return factors, {"q_matrix": t1 - t0, "kernel": t2 - t1, "split": t3 - t2}


METHODS: dict[str, Splitter] = {
    "berlekamp_random": lambda f, rng: berlekamp(f, rng, deterministic=False)[0],
    "berlekamp_deterministic": lambda f, rng: berlekamp(f, rng, deterministic=True)[0],
    "cantor_zassenhaus": cantor_zassenhaus,
    "sage_factor": lambda f, rng: [h for h, _ in f.factor()],
}


def random_squarefree(p: int, degree: int, rng: random.Random) -> Polynomial:
    ring = PolynomialRing(GF(p), "x")
    while True:
        f = ring([rng.randrange(p) for _ in range(degree)] + [1])
        if f.gcd(f.derivative()).degree() == 0:
            return f


def time_method(method: str, f: Polynomial, rng: random.Random) -> tuple[float, bool]:
    start = time.perf_counter()
    factors = METHODS[method](f, rng)
    elapsed = time.perf_counter() - start
    expected = sorted(tuple(h.list()) for h, _ in f.factor())
    return elapsed, sorted(tuple(h.list()) for h in factors) == expected


def measure(methods: list[str], p: int, degree: int, rng: random.Random) -> list[dict]:
    """Time every method on the same TRIALS random polynomials, so rows are directly comparable."""
    inputs = [random_squarefree(p, degree, rng) for _ in range(TRIALS)]
    mean_factors = statistics.mean(len(f.factor()) for f in inputs)
    rows = []
    for method in methods:
        results = [time_method(method, f, rng) for f in inputs]
        rows.append({
            "method": method, "p_bits": p.bit_length(), "p": p, "degree": degree,
            "median_s": statistics.median(t for t, _ in results), "trials": TRIALS,
            "verified": sum(ok for _, ok in results), "mean_factors": mean_factors,
        })
    return rows


def sweep_degree(rng: random.Random) -> list[dict]:
    methods = ["berlekamp_random", "cantor_zassenhaus", "sage_factor"]
    return [row for d in DEGREES for row in measure(methods, LARGE_PRIME, d, rng)]


def sweep_prime(rng: random.Random) -> list[dict]:
    rows = []
    for bits in RANDOMIZED_BITS:
        methods = ["berlekamp_random", "cantor_zassenhaus"]
        if bits in DETERMINISTIC_BITS:
            methods.insert(0, "berlekamp_deterministic")
        rows += measure(methods, int(previous_prime(2**bits)), FIXED_DEGREE, rng)
    return rows


def sweep_stages(rng: random.Random) -> list[dict]:
    rows = []
    for d in DEGREES:
        stages = [berlekamp(random_squarefree(LARGE_PRIME, d, rng), rng, deterministic=False)[1]
                  for _ in range(TRIALS)]
        rows.append({"degree": d, **{k: statistics.median(s[k] for s in stages) for k in STAGE_COLORS}})
    return rows


def write_csv(rows: list[dict], name: str) -> None:
    with open(RESULTS / name, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _style(ax, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, color=INK, loc="left", fontsize=13)
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color="#e4e3df", linewidth=0.8)
    ax.tick_params(colors=MUTED)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c9c8c2")


def plot_lines(rows: list[dict], x_key: str, title: str, xlabel: str, name: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8), facecolor=SURFACE)
    for method in COLORS:
        series = sorted((r[x_key], r["median_s"]) for r in rows if r["method"] == method)
        if not series:
            continue
        xs, ys = zip(*series)
        ax.plot(xs, ys, color=COLORS[method], linewidth=2, marker="o", markersize=5, label=LABELS[method])
        ax.annotate(LABELS[method], (xs[-1], ys[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=9, color=INK)
    ax.set_yscale("log")
    if x_key == "p_bits":
        ax.set_xscale("log", base=2)
        ticks = sorted({r["p_bits"] for r in rows if r["p_bits"] & (r["p_bits"] - 1) == 0})
        ax.set_xticks(ticks, [str(t) for t in ticks])
    _style(ax, title, xlabel, "median time per polynomial (s, log scale)")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="upper left")
    fig.tight_layout()
    fig.savefig(RESULTS / name, dpi=160)
    plt.close(fig)


def plot_stages(rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8), facecolor=SURFACE)
    names = {"q_matrix": "build Q", "kernel": "kernel of Q - I", "split": "randomized split"}
    bottom = [0.0] * len(rows)
    for stage, color in STAGE_COLORS.items():
        values = [r[stage] for r in rows]
        ax.bar([r["degree"] for r in rows], values, bottom=bottom, width=7, color=color,
               edgecolor=SURFACE, linewidth=2, label=names[stage])
        bottom = [b + v for b, v in zip(bottom, values)]
    _style(ax, f"Where Berlekamp spends its time (p = {LARGE_PRIME})", "degree d", "median time (s)")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="upper left")
    fig.tight_layout()
    fig.savefig(RESULTS / "berlekamp_stages.png", dpi=160)
    plt.close(fig)


def write_summary(by_degree: list[dict], by_prime: list[dict], stages: list[dict]) -> None:
    checked = [r for r in by_degree + by_prime if r["method"] != "sage_factor"]
    lines = [
        "# Benchmark results",
        "",
        f"Generated by `python -m bench.bench`. Each cell is the median over the same {TRIALS} random",
        "square-free monic polynomials for every method in that row. Every output of our code was",
        "checked against Sage `factor()`:",
        f"**{sum(r['verified'] for r in checked)} of {sum(r['trials'] for r in checked)} factorizations correct.**",
        "",
        f"## Time vs degree (p = {LARGE_PRIME})",
        "",
        "| d | " + " | ".join(LABELS[m] for m in ["berlekamp_random", "cantor_zassenhaus", "sage_factor"]) + " |",
        "|---|---|---|---|",
    ]
    for d in DEGREES:
        cells = [_seconds(r) for m in ["berlekamp_random", "cantor_zassenhaus", "sage_factor"]
                 for r in by_degree if r["degree"] == d and r["method"] == m]
        lines.append(f"| {d} | " + " | ".join(cells) + " |")
    lines += ["", f"## Time vs size of p (d = {FIXED_DEGREE})", "",
              "| bits of p | Berlekamp (deterministic) | Berlekamp (randomized) | Cantor-Zassenhaus (HW3) |",
              "|---|---|---|---|"]
    for b in sorted({r["p_bits"] for r in by_prime}):
        cells = []
        for m in ["berlekamp_deterministic", "berlekamp_random", "cantor_zassenhaus"]:
            match = [r for r in by_prime if r["p_bits"] == b and r["method"] == m]
            cells.append(_seconds(match[0]) if match else "-")
        lines.append(f"| {b} | " + " | ".join(cells) + " |")
    lines += ["", "Sage stores polynomials with FLINT for p < 2^63 and with NTL above, so part of the jump",
              "between 32 and 64 bits is the library, not the algorithm. Compare methods within a row.",
              "Deterministic Berlekamp is skipped above 16 bits: it tries all p values of s."]
    lines += ["", f"## Berlekamp stages (p = {LARGE_PRIME})", "",
              "| d | build Q | kernel of Q - I | randomized split |", "|---|---|---|---|"]
    lines += [f"| {r['degree']} | {r['q_matrix']:.4f} | {r['kernel']:.4f} | {r['split']:.4f} |" for r in stages]
    (RESULTS / "summary.md").write_text("\n".join(lines) + "\n")


def _seconds(row: dict) -> str:
    return f"{row['median_s']:.4f}"


def _quiet(_: str) -> None:
    pass


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    rng = random.Random(2026)
    by_degree = sweep_degree(rng)
    by_prime = sweep_prime(rng)
    stages = sweep_stages(rng)
    write_csv(by_degree, "time_vs_degree.csv")
    write_csv(by_prime, "time_vs_p.csv")
    write_csv(stages, "berlekamp_stages.csv")
    plot_lines(by_degree, "degree", f"Time vs degree (p = {LARGE_PRIME})", "degree d", "time_vs_degree.png")
    plot_lines(by_prime, "p_bits", f"Time vs size of p (d = {FIXED_DEGREE})", "bits of p (log scale)",
               "time_vs_p.png")
    plot_stages(stages)
    write_summary(by_degree, by_prime, stages)
    print((RESULTS / "summary.md").read_text())


if __name__ == "__main__":
    main()
