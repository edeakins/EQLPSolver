#!/usr/bin/env python3
"""
Prototype for the nonsingular combinatorial matching (companion to
matching_k4_findings.tex, section "Toward a fix").

Models the K4 fractional-matching instance's post-lift reduced basis and
implements:
  * ground truth   : exact rational rank of the reduced basis
  * Oracle A        : conservative forest/triangular test (provably safe, any LP)
  * Oracle B        : parity/even-cycle test (exact for uniform blocks in the
                      full-tight regime; has false-accepts on partial matchings)

Run:  python3 matching_oracle_prototype.py

Key results reproduced:
  - Oracle A: 0 false accepts over all 46 K4 candidate matchings, but places 0
    residuals on K4 (safely degrades to baseline).
  - Oracle B greedy: places {c2<-r2, c3<-r4, c4<-r6}, a rank-6 NONSINGULAR
    size-3 matching (3 pivots saved vs baseline), but is not unconditionally
    safe (see mismatch report).
"""
from fractions import Fraction as F
from itertools import product

# ---- K4 fractional-matching instance ---------------------------------------
# edges x1..x6 (0-based) ; residual r_i = x_i - x1 (parent x1), i = 2..6
inc = {'c1': {0, 1, 2}, 'c2': {0, 3, 4}, 'c3': {1, 3, 5}, 'c4': {2, 4, 5}}
resid_var = {2: 1, 3: 2, 4: 3, 5: 4, 6: 5}          # residual i -> edge index
adj = {'c2': [2, 3, 6], 'c3': [2, 4, 6], 'c4': [3, 5, 6]}  # slack ~ residual


def row_constraint(c):
    v = [F(0)] * 6
    for e in inc[c]:
        v[e] = F(1)
    return v


def row_link(i):                                     # x_i - x1
    v = [F(0)] * 6
    v[resid_var[i]] = F(1)
    v[0] -= F(1)
    return v


def reduced_basis(matching):
    """matching: list of (slack, residual). Returns the reduced structural rows
    (basic residual columns and their linking rows removed -- the A-submatrix)."""
    tight = ['c1'] + [s for s, _ in matching]
    cut = {r for _, r in matching}
    return [row_constraint(c) for c in tight] + \
           [row_link(i) for i in resid_var if i not in cut]


def rank(M):
    M = [r[:] for r in M]
    m, pr, rk = len(M), 0, 0
    for c in range(6):
        s = next((r for r in range(pr, m) if M[r][c] != 0), None)
        if s is None:
            continue
        M[pr], M[s] = M[s], M[pr]
        pv = M[pr][c]
        for r in range(m):
            if r != pr and M[r][c] != 0:
                f = M[r][c] / pv
                M[r] = [M[r][k] - f * M[pr][k] for k in range(6)]
        pr += 1
        rk += 1
    return rk


def truth_nonsingular(matching):
    R = reduced_basis(matching)
    return rank(R) == len(R)


# ---- structural transversal (bipartite matching on the sparsity pattern) ----
def has_transversal(M):
    n = len(M)
    if n != 6:
        return False
    A = [[c for c in range(6) if M[r][c] != 0] for r in range(n)]
    mc = [-1] * 6

    def aug(r, seen):
        for c in A[r]:
            if seen[c]:
                continue
            seen[c] = True
            if mc[c] == -1 or aug(mc[c], seen):
                mc[c] = r
                return True
        return False
    return sum(aug(r, [False] * 6) for r in range(6)) == 6


# ---- Oracle A: forest / triangular (provably safe) --------------------------
def is_forest(M):
    n = len(M)
    p = list(range(n + 6))

    def find(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    for r in range(n):
        for c in range(6):
            if M[r][c] != 0:
                a, b = find(r), find(n + c)
                if a == b:
                    return False          # cycle -> not a forest
                p[a] = b
    return True


def oracle_A(matching):
    M = reduced_basis(matching)
    return has_transversal(M) and is_forest(M)


# ---- Oracle B: parity union-find / even-cycle (uniform, full-tight) ---------
Gedge = {'x1': ('c1', 'c2'), 'x2': ('c1', 'c3'), 'x3': ('c1', 'c4'),
         'x4': ('c2', 'c3'), 'x5': ('c2', 'c4'), 'x6': ('c3', 'c4')}
var_of_resid = {i: f"x{resid_var[i] + 1}" for i in resid_var}


class ParityUF:
    def __init__(self, nodes):
        self.p = {n: n for n in nodes}
        self.par = {n: 0 for n in nodes}

    def find(self, x):
        if self.p[x] == x:
            return x, 0
        r, pr = self.find(self.p[x])
        self.p[x] = r
        self.par[x] ^= pr
        return r, self.par[x]

    def union(self, a, b, rel):
        ra, pa = self.find(a)
        rb, pb = self.find(b)
        if ra == rb:
            return (pa ^ pb) == rel
        self.p[ra] = rb
        self.par[ra] = pa ^ pb ^ rel
        return True


def oracle_B(matching):
    tight = ['c1'] + [s for s, _ in matching]
    T = set(tight)
    cutvars = {var_of_resid[r] for _, r in matching}
    if not has_transversal(reduced_basis(matching)):
        return False
    puf = ParityUF(tight)
    for v in cutvars:                      # cut edges must CROSS the signing
        a, b = Gedge[v]
        if a in T and b in T:
            if not puf.union(a, b, 1):     # odd cycle among cut edges -> safe
                return True
    comp, parity = {}, {}
    for c in tight:
        r, pp = puf.find(c)
        comp.setdefault(r, []).append(c)
        parity[c] = pp
    uncut = [v for v in Gedge if v not in cutvars]
    d = {c: 0 for c in tight}
    for v in uncut:
        a, b = Gedge[v]
        if a in T:
            d[a] += 1
        if b in T:
            d[b] += 1
    Ws = [sum((1 if parity[c] == 0 else -1) * d[c] for c in mem)
          for mem in comp.values()]
    for signs in product([1, -1], repeat=len(Ws)):
        if sum(s * w for s, w in zip(signs, Ws)) == 0:
            return False                   # even-cycle dependency
    return True


# ---- exhaustive validation + greedy demo -----------------------------------
def all_matchings():
    slacks = ['c2', 'c3', 'c4']
    seen, out = set(), []

    def rec(us, ur, cur):
        yield list(cur)
        for s in slacks:
            if s in us:
                continue
            for r in adj[s]:
                if r in ur:
                    continue
                yield from rec(us | {s}, ur | {r}, cur + [(s, r)])
    for m in rec(set(), set(), []):
        k = tuple(sorted(m))
        if k not in seen:
            seen.add(k)
            out.append(m)
    return out


def greedy(oracle, label):
    matching, used = [], set()
    print(f"\nGreedy with {label}:")
    for s in ['c2', 'c3', 'c4']:
        placed = False
        for r in adj[s]:
            if r in used:
                continue
            if oracle(matching + [(s, r)]):
                matching.append((s, r))
                used.add(r)
                print(f"  accept (c{s[1]} <- r{r})")
                placed = True
                break
            else:
                print(f"  reject (c{s[1]} <- r{r})")
        if not placed:
            print(f"  {s}: nothing placeable")
    print(f"  FINAL {matching} (size {len(matching)}), "
          f"truth nonsingular: {truth_nonsingular(matching)}")


if __name__ == "__main__":
    ms = all_matchings()
    faA = sum(1 for m in ms if oracle_A(m) and not truth_nonsingular(m))
    faB = sum(1 for m in ms if oracle_B(m) and not truth_nonsingular(m))
    print(f"candidate matchings: {len(ms)}")
    print(f"Oracle A false-accepts (unsafe): {faA}   <- must be 0")
    print(f"Oracle B false-accepts (unsafe): {faB}   <- nonzero: not "
          f"unconditionally safe")
    greedy(oracle_A, "Oracle A (conservative forest)")
    greedy(oracle_B, "Oracle B (parity / even-cycle)")
