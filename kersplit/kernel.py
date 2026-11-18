from sage.matrix.matrix0 import Matrix


def left_kernel_basis(m: Matrix) -> list[list]:
    """Return a basis of {v : v m = 0}, by Gauss-Jordan elimination on the transpose of m."""
    field = m.base_ring()
    rows = [list(column) for column in m.columns()]
    n = m.nrows()
    pivot_columns = _reduce_to_rref(rows, n)
    basis = []
    for free in (c for c in range(n) if c not in pivot_columns):
        v = [field.zero()] * n
        v[free] = field.one()
        for row, pivot in zip(rows, pivot_columns):
            v[pivot] = -row[free]
        basis.append(v)
    return basis


def _reduce_to_rref(rows: list[list], n_columns: int) -> list[int]:
    """Reduce rows in place to reduced row echelon form and return the pivot columns."""
    pivot_columns: list[int] = []
    for col in range(n_columns):
        r = len(pivot_columns)
        pivot = next((i for i in range(r, len(rows)) if rows[i][col] != 0), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inverse = 1 / rows[r][col]
        rows[r] = [a * inverse for a in rows[r]]
        for i, row in enumerate(rows):
            if i != r and row[col] != 0:
                k = row[col]
                rows[i] = [a - k * b for a, b in zip(row, rows[r])]
        pivot_columns.append(col)
    return pivot_columns
