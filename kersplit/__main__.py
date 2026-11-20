import argparse
import sys

from kersplit.berlekamp import METHODS, factor
from kersplit.polyio import parse_problem


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="python -m kersplit",
        description="Factor f in Z_p[x] with Berlekamp's algorithm. Reads HW3-format input on stdin.",
    )
    parser.add_argument("--method", choices=METHODS, default="auto")
    parser.add_argument("--trace", action="store_true", help="print every step")
    parser.add_argument("--seed", type=int, default=0, help="seed for the random split")
    args = parser.parse_args()
    try:
        problem = parse_problem(sys.stdin.read())
        result = factor(problem.f, method=args.method, trace=args.trace, seed=args.seed)
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    if not args.trace:
        print(result)
    if len(result.factors) == 1 and result.factors[0][1] == 1:
        print("f is irreducible")
    return 0


if __name__ == "__main__":
    sys.exit(main())
