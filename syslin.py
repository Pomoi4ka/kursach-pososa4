def rref_with_pivots(A, tol=1e-12):
    m = len(A)
    n = len(A[0]) if m else 0
    # make copy
    M = [row[:] for row in A]
    row = 0
    pivots = []
    for col in range(n):
        # find pivot row
        sel = max(range(row, m), key=lambda i: abs(M[i][col]))
        if abs(M[sel][col]) <= tol:
            continue
        # swap
        M[row], M[sel] = M[sel], M[row]
        # normalize pivot row
        piv = M[row][col]
        M[row] = [val / piv for val in M[row]]
        print("M[row]: ", M[row])
        # eliminate other rows
        for i in range(m):
            if i == row: continue
            factor = M[i][col]
            if factor:
                M[i] = [M[i][j] - factor * M[row][j] for j in range(n)]
        pivots.append((row, col))
        row += 1
        print(M)
        if row == m: break
    return M, pivots

def null_space(A, tol=1e-12):
    m = len(A)
    n = len(A[0]) if m else 0
    # get RREF and pivot positions
    R, pivots = rref_with_pivots(A, tol=tol)
    pivot_cols = {col for _, col in pivots}
    free_cols = [j for j in range(n) if j not in pivot_cols]
    if not free_cols:
        return []  # only trivial solution
    basis = []
    for free in free_cols:
        vec = [0.0]*n
        vec[free] = 1.0
        # set pivot variables according to RREF: for each pivot row r with pivot col c,
        # x_c = - sum_{free j} R[r][j]*x_j  (here only x_free = 1)
        for r, c in pivots:
            vec[c] = -R[r][free]
        basis.append(vec)
    print(f"R = {R}")
    return basis

def null_space_basis(A, tol=1e-12):
    # (Use the null_space function from earlier assistant message.)
    # Returns list of basis vectors (each a list of length n).
    return null_space(A, tol=tol)

# Build a particular solution from coefficients c (list length k)
def build_solution(basis, c):
    n = len(basis[0]) if basis else 0
    x = [0.0]*n
    for coeff, vec in zip(c, basis):
        for i in range(n):
            x[i] += coeff * vec[i]
    return x

# Example: take coefficients c1=2, c2=-1

# Example usage
if __name__ == "__main__":
    A = [
#        [1.0, 2.0, 8.0, 1.4],
#        [1.0, 2.5, 1.0, 1],
#        [6.0, 4.0, 3.0, 1],
        [1, 2],
        [1, 5],
    ]
    print(*A, sep='\n')
    basis = null_space(A)
    print("Null-space basis (each vector is a column):")
    for v in basis:
        print(v)
    if basis:
        x = build_solution(basis, [1.0]*len(basis))  # example coefficients
        print("Example solution:", x)
    else:
        print("Only trivial solution x = 0")
    # Expected: two basis vectors spanning the nullspace (since rank=2, n=4 -> nullity=2)
