#!/usr/bin/env python3
"""
gen_soft_edge_patterns_p2.py

Soft-decision (Chase p=2) EDGE-CASE pattern generator for BCH modes:
  code=1: (63,51)   m=6  t=2   pattern IDs start at 450
  code=2: (255,239) m=8  t=2   pattern IDs start at 550
  code=3: (1023,983)m=10 t=4   pattern IDs start at 650

This generates "most difficult" inputs for soft decoding, stressing:
  - Errors in a row (bursts)
  - Errors at front/end (0, n-1, n-2,...)
  - Least reliable bits forced to be at tail positions (n-2, n-1) or head positions
  - v in [0, t+p] where p=2; if v>t, extra errors are among the 2 least reliable bits
  - LLR in [-127,127], no -128
  - One codeword per pattern
  - 4 Chase variants per pattern (flip/not flip 2 least reliable bits)

Decoding correctness:
  - Uses normalized BM to compute sigma/roots (reference).
  - HR shows "sigma_raw" as a scalar multiple of sigma_norm (representation-only).
  - Answer file outputs FINAL correction locations relative to original hard word:
        e_total = flipmask XOR roots
    (so if you flipped LR bits that were errors, they appear in the output.)

Outputs per pattern:
  pXXX.txt   : packed LLR stream
  pXXXa.txt  : winner e_total locations (10-bit unsigned bin), or 1023 if empty
  pXXXhr.txt : details incl LLRs, true errors, S1..S2t, sigma_norm/sigma_raw,
               roots, flip_locs, e_total, decoding_failed, correlation.

"""

from __future__ import annotations
import argparse
import os
import random
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional


# -------------------------
# GF(2^m) arithmetic
# -------------------------
@dataclass(frozen=True)
class GF2m:
    m: int
    prim_poly: int
    exp: List[int]
    log: List[int]

    @staticmethod
    def build(m: int, prim_poly: int) -> "GF2m":
        size = 1 << m
        exp = [0] * (size - 1)
        log = [-1] * size

        def mul(a: int, b: int) -> int:
            res = 0
            aa = a
            bb = b
            mask = (1 << m) - 1
            red = prim_poly & mask
            while bb:
                if bb & 1:
                    res ^= aa
                bb >>= 1
                carry = aa & (1 << (m - 1))
                aa = (aa << 1) & mask
                if carry:
                    aa ^= red
            return res & mask

        exp[0] = 1
        for k in range(1, size - 1):
            exp[k] = mul(exp[k - 1], 2)
        for k, v in enumerate(exp):
            log[v] = k
        return GF2m(m=m, prim_poly=prim_poly, exp=exp, log=log)

    def mul(self, a: int, b: int) -> int:
        res = 0
        aa = a
        bb = b
        mask = (1 << self.m) - 1
        red = self.prim_poly & mask
        while bb:
            if bb & 1:
                res ^= aa
            bb >>= 1
            carry = aa & (1 << (self.m - 1))
            aa = (aa << 1) & mask
            if carry:
                aa ^= red
        return res & mask

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("GF inverse of 0")
        N = len(self.exp)
        return self.exp[(N - self.log[a]) % N]

    def alpha_pow(self, k: int) -> int:
        return self.exp[k % len(self.exp)]

    def to_alpha_str(self, a: int) -> str:
        if a == 0:
            return "0"
        if a == 1:
            return "1"
        return f"a^{self.log[a]}"

    def to_bits(self, a: int) -> str:
        return format(a, f"0{self.m}b")


# -------------------------
# BCH mode config
# -------------------------
@dataclass(frozen=True)
class BCHMode:
    code: int
    n: int
    k: int
    t: int
    m: int
    prim_poly: int

    @property
    def llr_count(self) -> int:
        return self.n + 1

    @property
    def lines_per_codeword(self) -> int:
        assert self.llr_count % 8 == 0
        return self.llr_count // 8


MODES = {
    1: BCHMode(code=1, n=63,   k=51,  t=2, m=6,  prim_poly=0x43),
    2: BCHMode(code=2, n=255,  k=239, t=2, m=8,  prim_poly=0x11D),
    3: BCHMode(code=3, n=1023, k=983, t=4, m=10, prim_poly=0x409),
}


# -------------------------
# IO formatting helpers
# -------------------------
def int8_to_twos_comp_bits(x: int) -> str:
    if x == -128:
        raise ValueError("Soft decision forbids -128")
    if not (-127 <= x <= 127):
        raise ValueError(f"int8 out of range [-127,127]: {x}")
    if x < 0:
        x = (1 << 8) + x
    return format(x, "08b")


def pack_8_llrs_to_line(llrs: List[int]) -> str:
    assert len(llrs) == 8
    return "".join(int8_to_twos_comp_bits(v) for v in llrs)


def u10_bits(x: int) -> str:
    if not (0 <= x <= 1023):
        raise ValueError(f"u10 out of range: {x}")
    return format(x, "010b")


def llr_stream_from_l_i(l_i: List[int]) -> List[int]:
    """Stream order LLR0..LLRn with LLR1->r_{n-1}, ..., LLRn->r_0."""
    n = len(l_i)
    llrs = [0] * (n + 1)
    llrs[0] = 0
    for sidx in range(1, n + 1):
        i = n - sidx
        llrs[sidx] = l_i[i]
    return llrs


def hard_bits_from_l_i(l_i: List[int]) -> List[int]:
    return [0 if li >= 0 else 1 for li in l_i]


# -------------------------
# BCH core: syndromes, BM, Chien
# -------------------------
def syndromes(gf: GF2m, r_bits: List[int], t: int) -> List[int]:
    n = len(r_bits)
    S = []
    for k in range(1, 2 * t + 1):
        acc = 0
        for i, b in enumerate(r_bits):
            if b:
                acc ^= gf.alpha_pow(k * i)
        S.append(acc)
    return S


def berlekamp_massey_normalized(gf: GF2m, S: List[int], t: int) -> List[int]:
    C = [0] * (t + 1)
    B = [0] * (t + 1)
    C[0] = 1
    B[0] = 1
    L = 0
    m = 1
    b = 1

    for n_idx in range(0, 2 * t):
        d = S[n_idx]
        for i in range(1, L + 1):
            if C[i] and S[n_idx - i]:
                d ^= gf.mul(C[i], S[n_idx - i])

        if d == 0:
            m += 1
            continue

        T = C[:]
        coef = gf.mul(d, gf.inv(b))
        for i in range(0, t + 1 - m):
            if B[i]:
                C[i + m] ^= gf.mul(coef, B[i])

        if 2 * L <= n_idx:
            B = T
            L = n_idx + 1 - L
            b = d
            m = 1
        else:
            m += 1

    sigma = C[:L + 1] if L > 0 else [1]
    if sigma[0] != 1 and sigma[0] != 0:
        inv0 = gf.inv(sigma[0])
        sigma = [gf.mul(c, inv0) for c in sigma]
    elif sigma[0] == 0:
        sigma = [1]
    return sigma


def scale_sigma_as_raw(gf: GF2m, sigma_norm: List[int], rng: random.Random) -> Tuple[int, List[int]]:
    k = rng.randint(1, (1 << gf.m) - 1)
    return k, [gf.mul(k, c) for c in sigma_norm]


def poly_eval(gf: GF2m, poly: List[int], x: int) -> int:
    acc = 0
    xp = 1
    for c in poly:
        if c:
            acc ^= gf.mul(c, xp)
        xp = gf.mul(xp, x)
    return acc


def chien_search(gf: GF2m, sigma: List[int], n: int,
                 root_at_inverse: bool = True,
                 loc_is_degree: bool = True) -> List[int]:
    roots = []
    for j in range(n):
        x = gf.alpha_pow((-j) % n) if root_at_inverse else gf.alpha_pow(j)
        if poly_eval(gf, sigma, x) == 0:
            loc = j if loc_is_degree else (n - 1 - j)
            roots.append(loc)
    roots.sort()
    return roots


def correlation(l_i: List[int], c_hat: List[int]) -> int:
    s = 0
    for li, ci in zip(l_i, c_hat):
        s += li * (1 if ci == 0 else -1)
    return s


def xor_locations(a: List[int], b: List[int]) -> List[int]:
    s = set(a)
    for x in b:
        if x in s:
            s.remove(x)
        else:
            s.add(x)
    return sorted(s)


# -------------------------
# Soft EDGE case construction
# -------------------------
@dataclass(frozen=True)
class EdgeCaseSpec:
    name: str
    # desired least-reliable positions (2 positions)
    lr_pos: Tuple[int, int]
    # total errors v (0..t+p)
    v: int
    # how many of the errors are among lr_pos (must be max(0,v-t))
    extra_in_lr: int
    # other error locations (besides lr errors), may include bursts/ends etc.
    other_errors: List[int]
    # whether to enforce a burst around a pivot (for auto-fill)
    burst_len: int = 0
    burst_start: Optional[int] = None


def unique_sorted_in_range(xs: List[int], n: int) -> List[int]:
    return sorted({x for x in xs if 0 <= x < n})


def make_edge_case(mode: BCHMode, spec: EdgeCaseSpec, rng: random.Random, p: int = 2) -> Tuple[List[int], List[int], Tuple[int, int], str]:
    """
    Produce (l_i, e_true, lr_pos, tag) for the given edge spec.
    Enforces:
      - v in [0, t+p]
      - if v>t, extra errors are within lr_pos (extra_in_lr must match)
      - lr_pos are the two smallest |LLR| (magnitudes 1 and 2)
      - all other magnitudes >=10
    """
    n, t = mode.n, mode.t
    assert p == 2
    assert 0 <= spec.v <= t + p
    assert 0 <= spec.extra_in_lr <= 2

    lr0, lr1 = spec.lr_pos
    if lr0 == lr1 or not (0 <= lr0 < n and 0 <= lr1 < n):
        raise ValueError(f"Bad lr_pos {spec.lr_pos} for n={n}")
    lr_pos = (lr0, lr1) if lr0 < lr1 else (lr1, lr0)

    # Determine how many errors must be in lr_pos
    required_extra = max(0, spec.v - t)
    if spec.extra_in_lr != required_extra:
        raise ValueError(f"{spec.name}: extra_in_lr must equal max(0,v-t)={required_extra}")

    e_set = set()

    # Put required extra errors into LR positions:
    if required_extra == 1:
        # choose one (prefer lr0 to be deterministic if it is in other_errors, else random)
        chosen = lr_pos[0] if lr_pos[0] in spec.other_errors else rng.choice(list(lr_pos))
        e_set.add(chosen)
    elif required_extra == 2:
        e_set.add(lr_pos[0])
        e_set.add(lr_pos[1])

    # Add other specified errors
    e_set |= set(spec.other_errors)

    # Optional burst fill
    if spec.burst_len and spec.burst_start is not None:
        for k in range(spec.burst_len):
            e_set.add(spec.burst_start + k)

    # Trim/extend to exactly v errors
    e_list = unique_sorted_in_range(list(e_set), n)

    # If too many, trim but keep LR-required ones and endpoints/burst preference
    if len(e_list) > spec.v:
        must = set()
        # keep required LR errors
        if required_extra == 1:
            must |= (set(e_list) & set(lr_pos))
        elif required_extra == 2:
            must |= set(lr_pos)
        # keep hard endpoints if present
        for x in [0, 1, n - 2, n - 1]:
            if x in e_list:
                must.add(x)
        # fill remaining from e_list order
        kept = []
        for x in e_list:
            if x in must and x not in kept:
                kept.append(x)
        for x in e_list:
            if len(kept) >= spec.v:
                break
            if x not in kept:
                kept.append(x)
        e_list = sorted(kept)

    # If too few, add more near edges (harder) excluding LR if already used
    if len(e_list) < spec.v:
        candidates = [0, 1, 2, 3, n - 1, n - 2, n - 3, n - 4, n // 2, n // 2 - 1, n // 2 + 1]
        for c in candidates:
            if 0 <= c < n and c not in e_list:
                e_list.append(c)
                if len(e_list) == spec.v:
                    break
        # still too few: random fill
        if len(e_list) < spec.v:
            pool = [i for i in range(n) if i not in e_list]
            rng.shuffle(pool)
            while len(e_list) < spec.v and pool:
                e_list.append(pool.pop())
        e_list = sorted(e_list)

    e_true = e_list
    assert len(e_true) == spec.v

    # Build r bits (tx all-zero)
    r_bits = [0] * n
    for j in e_true:
        r_bits[j] = 1

    # Build LLRs:
    #  - lr positions get magnitudes 1 and 2 (smallest |LLR|).
    #  - optionally make lr positions sit at the last inputs => spec.lr_pos handles that.
    l_i = [0] * n
    mags_lr = [1, 2]
    for idx, pos in enumerate(lr_pos):
        mag = mags_lr[idx]
        l_i[pos] = -mag if r_bits[pos] else +mag

    # other positions: make "hard" ones huge, but also make some tricky near-threshold
    # We'll mostly keep abs>=10; but for edge hardness, we sprinkle a few abs=10..15
    # near the error burst region to stress correlation.
    for i in range(n):
        if i in lr_pos:
            continue
        if i in e_true:
            # error bits: negative; choose magnitude moderate (10..40) to avoid being super obvious
            mag = rng.randint(10, 40)
            l_i[i] = -mag
        else:
            # non-error bits: positive; some are weak (10..20), others strong (60..127)
            if rng.random() < 0.15:
                mag = rng.randint(10, 20)
            else:
                mag = rng.randint(60, 127)
            l_i[i] = +mag

    tag = spec.name
    return l_i, e_true, lr_pos, tag


def edge_suite_for_mode(mode: BCHMode) -> List[EdgeCaseSpec]:
    """
    Return a deterministic list of edge-case specs for a given mode.
    IDs will cycle through these specs as you generate multiple patterns.
    """
    n, t = mode.n, mode.t
    p = 2

    # Helpers for LR positions at ends
    LR_TAIL = (n - 2, n - 1)
    LR_HEAD = (0, 1)
    LR_MIX  = (0, n - 1)

    suite: List[EdgeCaseSpec] = []

    # Case A: v=t+p, extra errors are BOTH LR at tail, remaining errors are a burst at start
    # For t=2: v=4 -> errors: {tail LR} + {0,1} burst
    # For t=4: v=6 -> errors: {tail LR} + {0,1,2,3} burst
    suite.append(EdgeCaseSpec(
        name="v=max, LR=tail, burst@head",
        lr_pos=LR_TAIL,
        v=t + p,
        extra_in_lr=2,  # since v>t always for v=t+2
        other_errors=list(range(0, t))  # burst length t at head
    ))

    # Case B: v=t+p, LR at head, remaining errors a burst at tail (excluding LR if overlaps)
    other_tail = list(range(max(0, n - t), n))
    suite.append(EdgeCaseSpec(
        name="v=max, LR=head, burst@tail",
        lr_pos=LR_HEAD,
        v=t + p,
        extra_in_lr=2,
        other_errors=other_tail
    ))

    # Case C: v=t+p, LR mixed at (0,n-1), remaining errors in middle consecutive burst
    mid = n // 2
    burst_mid = list(range(mid - (t // 2), mid - (t // 2) + t))
    suite.append(EdgeCaseSpec(
        name="v=max, LR=(0,n-1), burst@mid",
        lr_pos=LR_MIX,
        v=t + p,
        extra_in_lr=2,
        other_errors=burst_mid
    ))

    # Case D: v=t+1 (only one extra), LR at tail, extra error is exactly n-1, other errors are consecutive near head
    if t >= 2:
        suite.append(EdgeCaseSpec(
            name="v=t+1, extra in LR tail, burst@head",
            lr_pos=LR_TAIL,
            v=t + 1,
            extra_in_lr=1,
            other_errors=list(range(0, t))  # t errors at head; will be trimmed to v
        ))

    # Case E: v=t, but LR bits at tail are NOT errors (tests flipping wrong bits)
    # Here v==t => required_extra=0, but LR bits are very small magnitude and *correct*.
    # We choose errors as a burst at head or tail.
    suite.append(EdgeCaseSpec(
        name="v=t, LR=tail (NOT errors), errors=burst@head",
        lr_pos=LR_TAIL,
        v=t,
        extra_in_lr=0,
        other_errors=list(range(0, t))
    ))

    # Case F: v=t, errors at both ends + clustered (stress ends)
    ends = [0, n - 1]
    more = [1, n - 2]
    needed = ends + more
    suite.append(EdgeCaseSpec(
        name="v=t, errors near both ends",
        lr_pos=LR_TAIL,  # LR at tail, but v=t => LR need not be errors
        v=t,
        extra_in_lr=0,
        other_errors=needed[:t]
    ))

    # Case G: v=0, but LR bits at tail (tie-like but still valid)
    suite.append(EdgeCaseSpec(
        name="v=0, LR=tail",
        lr_pos=LR_TAIL,
        v=0,
        extra_in_lr=0,
        other_errors=[]
    ))

    # Case H: v=1, error at n-1, LR bits at tail (one LR is error)
    suite.append(EdgeCaseSpec(
        name="v=1, error at tail, LR=tail",
        lr_pos=LR_TAIL,
        v=1,
        extra_in_lr=0,  # v<=t, so no requirement; but we *choose* to place error at LR
        other_errors=[n - 1]
    ))

    # Case I: v=t+p, LR at tail, all errors consecutive at tail (burst ends include LR)
    tail_burst = list(range(n - (t + p), n))
    suite.append(EdgeCaseSpec(
        name="v=max, all errors consecutive at tail",
        lr_pos=LR_TAIL,
        v=t + p,
        extra_in_lr=2,
        other_errors=tail_burst
    ))

    return suite


# -------------------------
# Run Chase on one codeword, write pXXX files
# -------------------------
def write_one_soft_edge_pattern(outdir: str,
                                pid: int,
                                mode: BCHMode,
                                seed: int,
                                spec: EdgeCaseSpec,
                                p: int = 2,
                                root_at_inverse: bool = True,
                                loc_is_degree: bool = True) -> None:
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    gf = GF2m.build(mode.m, mode.prim_poly)

    l_i, e_true, lr_pos, tag = make_edge_case(mode, spec, rng, p=p)

    r_base = hard_bits_from_l_i(l_i)
    abs_sorted = sorted(range(mode.n), key=lambda i: (abs(l_i[i]), i))
    lrA, lrB = abs_sorted[0], abs_sorted[1]

    flip_masks = [(0, 0), (0, 1), (1, 0), (1, 1)]
    variants = []

    for vidx, (fa, fb) in enumerate(flip_masks):
        r_cand = r_base[:]
        flip_locs = []
        if fa:
            r_cand[lrA] ^= 1
            flip_locs.append(lrA)
        if fb:
            r_cand[lrB] ^= 1
            flip_locs.append(lrB)
        flip_locs.sort()

        S = syndromes(gf, r_cand, mode.t)
        sigma_norm = berlekamp_massey_normalized(gf, S, mode.t)
        deg = len(sigma_norm) - 1
        k_raw, sigma_raw = scale_sigma_as_raw(gf, sigma_norm, rng)

        roots = chien_search(
            gf, sigma_norm, mode.n,
            root_at_inverse=root_at_inverse,
            loc_is_degree=loc_is_degree
        ) if deg > 0 else []

        failed = (len(roots) < deg)

        if not failed:
            e_total_locs = xor_locations(flip_locs, roots)
            # decoded codeword for correlation: c_hat = r_base XOR e_total
            c_hat = r_base[:]
            for j in e_total_locs:
                c_hat[j] ^= 1
            corr = correlation(l_i, c_hat)
        else:
            e_total_locs = []
            corr = None

        variants.append({
            "vidx": vidx,
            "flip": (fa, fb),
            "flip_locs": flip_locs,
            "S": S,
            "sigma_norm": sigma_norm,
            "sigma_raw": sigma_raw,
            "k_raw": k_raw,
            "deg": deg,
            "roots": roots,
            "e_total_locs": e_total_locs,
            "failed": failed,
            "corr": corr,
        })

    successful = [v for v in variants if (not v["failed"]) and (v["corr"] is not None)]
    if successful:
        winner = max(successful, key=lambda v: v["corr"])
        win_total = winner["e_total_locs"]
    else:
        winner = None
        win_total = []

    fin = os.path.join(outdir, f"p{pid}.txt")
    fans = os.path.join(outdir, f"p{pid}a.txt")
    fhr  = os.path.join(outdir, f"p{pid}hr.txt")

    # input stream
    llr_stream = llr_stream_from_l_i(l_i)
    with open(fin, "w") as f_in:
        for li in range(mode.lines_per_codeword):
            chunk = llr_stream[li * 8:(li + 1) * 8]
            f_in.write(pack_8_llrs_to_line(chunk) + "\n")

    # answer: FINAL correction locations
    with open(fans, "w") as f_ans:
        if len(win_total) == 0:
            f_ans.write(u10_bits(1023) + "\n")
        else:
            for loc in win_total:
                f_ans.write(u10_bits(loc) + "\n")

    # HR
    with open(fhr, "w") as f:
        f.write(f"Pattern p{pid} (SOFT EDGE, p=2)\n")
        f.write(f"Mode code={mode.code}, (n,k)=({mode.n},{mode.k}), m={mode.m}, t={mode.t}\n")
        f.write(f"Primitive poly = 0x{mode.prim_poly:X}\n")
        f.write(f"Edge tag: {tag}\n")
        f.write(f"Chien convention: root_at_inverse={root_at_inverse}, loc_is_degree={loc_is_degree}\n\n")

        f.write("1) Original channel LLRs l_i for r_i (i=0..n-1, ascending degree):\n")
        f.write("   i : l_i (hard_bit)\n")
        for i in range(mode.n):
            hb = 0 if l_i[i] >= 0 else 1
            f.write(f"   {i:4d}: {l_i[i]:4d} ({hb})\n")
        f.write("\n")
        f.write(f"True error locations (tx all-zero): {e_true}\n\n")

        f.write("2) Two least reliable indices (by |LLR|):\n")
        f.write(f"   constructed lr_pos = {lr_pos}\n")
        f.write(f"   found by |LLR| sort = ({lrA},{lrB})\n")
        f.write(f"   lrA={lrA} l={l_i[lrA]} |l|={abs(l_i[lrA])}\n")
        f.write(f"   lrB={lrB} l={l_i[lrB]} |l|={abs(l_i[lrB])}\n\n")

        f.write("3/4/5) Chase variants:\n")
        f.write("   roots = e_hat locations relative to r_cand (after flips)\n")
        f.write("   e_total_locs = flipmask XOR roots (FINAL corrections relative to original r_base)\n")
        f.write("   decoding_failed = TRUE if (#roots < degree(sigma)).\n\n")

        for v in variants:
            mark = "  <== WINNER" if (winner is not None and winner["vidx"] == v["vidx"]) else ""
            fa, fb = v["flip"]
            f.write(f"--- Variant {v['vidx']} flip(lrA,lrB)=({fa},{fb}){mark}\n")
            f.write("    Syndromes S1..S2t:\n")
            for si, sval in enumerate(v["S"], start=1):
                f.write(f"      S{si}: bits={gf.to_bits(sval)}  alpha={gf.to_alpha_str(sval)}  int={sval}\n")

            f.write(f"    sigma degree (L) = {v['deg']}\n")
            f.write("    sigma_norm (BM, normalized):\n")
            for d, c in enumerate(v["sigma_norm"]):
                f.write(f"      sigma_norm[{d}]: bits={gf.to_bits(c)}  alpha={gf.to_alpha_str(c)}  int={c}\n")

            k = v["k_raw"]
            f.write(f"    sigma_raw (unnormalized display) = k*sigma_norm, k=bits={gf.to_bits(k)} alpha={gf.to_alpha_str(k)} int={k}\n")
            for d, c in enumerate(v["sigma_raw"]):
                f.write(f"      sigma_raw[{d}]: bits={gf.to_bits(c)}  alpha={gf.to_alpha_str(c)}  int={c}\n")

            f.write(f"    flip_locs = {v['flip_locs']}\n")
            f.write(f"    roots(e_hat) = {v['roots']}\n")
            f.write(f"    e_total_locs (flip XOR roots) = {v['e_total_locs']}\n")
            f.write(f"    decoding_failed = {v['failed']}\n")
            if v["corr"] is None:
                f.write("    correlation = N/A (failed)\n\n")
            else:
                f.write(f"    correlation = {v['corr']}\n\n")

        if winner is None:
            f.write("No successful variant found => answer outputs 1023.\n")
        else:
            f.write(f"Chosen winner variant = {winner['vidx']} with correlation = {winner['corr']}\n")
            f.write(f"Winner e_total_locs (FINAL output) = {winner['e_total_locs']}\n")

    print(f"[OK] wrote {fin}, {fans}, {fhr}")


# -------------------------
# Main
# -------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="patterns_soft_edge_p2", help="output directory")
    ap.add_argument("--num_patterns", type=int, default=10, help="patterns per mode (default: 10)")
    ap.add_argument("--seed", type=int, default=20251213, help="global seed")
    ap.add_argument("--chien_root_at_alpha_j", action="store_true", default=False,
                    help="Use sigma(alpha^j)==0 convention (default: sigma(alpha^{-j})==0)")
    ap.add_argument("--loc_is_msb", action="store_true", default=False,
                    help="Return location (n-1-j) instead of j")
    args = ap.parse_args()

    root_at_inverse = not args.chien_root_at_alpha_j
    loc_is_degree = not args.loc_is_msb

    base = {1: 450, 2: 550, 3: 650}

    for code, mode in MODES.items():
        suite = edge_suite_for_mode(mode)
        for pidx in range(args.num_patterns):
            pid = base[code] + pidx
            spec = suite[pidx % len(suite)]
            patt_seed = (args.seed * 1000003) ^ (pid * 9176) ^ (code * 1315423911)
            write_one_soft_edge_pattern(
                outdir=args.outdir,
                pid=pid,
                mode=mode,
                seed=patt_seed,
                spec=spec,
                p=2,
                root_at_inverse=root_at_inverse,
                loc_is_degree=loc_is_degree
            )


if __name__ == "__main__":
    main()
