#!/usr/bin/env python3
"""
gen_soft_patterns_p2_corrected_with_flipmask_in_answer.py

Soft-decision BCH Chase-pattern generator (p=2), corrected so that:
  - Decoding and roots are computed using normalized BM (reference-correct).
  - HR file shows sigma_raw as a scalar multiple of sigma_norm (representation-only).
  - Answer file outputs FINAL correction locations w.r.t. ORIGINAL hard word:
        e_total = e_flipmask XOR e_hat
    i.e., include flipped least-reliable indices in the output (as needed).

LLR restriction: [-127,127] (NO -128)
One codeword per pattern.

Pattern numbering:
  m=6  starts at 400
  m=8  starts at 500
  m=10 starts at 600
"""

from __future__ import annotations
import argparse
import os
import random
from dataclasses import dataclass
from typing import List, Tuple


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
# Packing helpers
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


# -------------------------
# Mapping: l_i (r_i) -> stream LLR[0..n]
# -------------------------
def llr_stream_from_l_i(l_i: List[int]) -> List[int]:
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
# Syndromes S1..S2t
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


# -------------------------
# Standard BM (normalized sigma)
# -------------------------
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
    sigma_raw = [gf.mul(k, c) for c in sigma_norm]
    return k, sigma_raw


# -------------------------
# Poly eval + Chien
# -------------------------
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
    for j in range(0, n):
        x = gf.alpha_pow((-j) % n) if root_at_inverse else gf.alpha_pow(j)
        if poly_eval(gf, sigma, x) == 0:
            loc = j if loc_is_degree else (n - 1 - j)
            roots.append(loc)
    roots.sort()
    return roots


# -------------------------
# Correlation
# -------------------------
def correlation(l_i: List[int], c_hat: List[int]) -> int:
    s = 0
    for li, ci in zip(l_i, c_hat):
        s += li * (1 if ci == 0 else -1)
    return s


# -------------------------
# Soft-case construction (p=2 constraints)
# -------------------------
def pick_two_least_reliable_positions(rng: random.Random, n: int) -> Tuple[int, int]:
    a = rng.randrange(n)
    b = rng.randrange(n - 1)
    if b >= a:
        b += 1
    return (a, b) if a < b else (b, a)


def make_soft_case(rng: random.Random, mode: BCHMode, p: int = 2) -> Tuple[List[int], List[int], Tuple[int, int]]:
    n, t = mode.n, mode.t
    v = rng.randint(0, t + p)

    lr0, lr1 = pick_two_least_reliable_positions(rng, n)
    lr_pos = (lr0, lr1)

    extra = max(0, v - t)  # 0..2
    e_set = set()

    if extra == 1:
        e_set.add(lr_pos[rng.randrange(2)])
    elif extra == 2:
        e_set.add(lr0)
        e_set.add(lr1)

    remain = v - len(e_set)
    if remain > 0:
        candidates = [i for i in range(n) if i not in e_set]
        for loc in rng.sample(candidates, remain):
            e_set.add(loc)

    e_true = sorted(e_set)
    r_bits = [0] * n
    for j in e_true:
        r_bits[j] = 1

    # LLRs: lr bits are abs=1,2; others abs>=10
    l_i = [0] * n
    mags_lr = [1, 2]
    for idx, pos in enumerate(lr_pos):
        mag = mags_lr[idx]
        l_i[pos] = -mag if r_bits[pos] == 1 else +mag

    for i in range(n):
        if i in lr_pos:
            continue
        mag = rng.randint(10, 127)
        l_i[i] = -mag if r_bits[i] == 1 else +mag

    return l_i, e_true, lr_pos


# -------------------------
# Pattern writer (one codeword)
# -------------------------
def xor_locations(a: List[int], b: List[int]) -> List[int]:
    """Return sorted locations for XOR of two location-sets."""
    s = set(a)
    for x in b:
        if x in s:
            s.remove(x)
        else:
            s.add(x)
    return sorted(s)


def write_one_soft_pattern(outdir: str,
                           pid: int,
                           mode: BCHMode,
                           seed: int,
                           p: int = 2,
                           root_at_inverse: bool = True,
                           loc_is_degree: bool = True) -> None:
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    gf = GF2m.build(mode.m, mode.prim_poly)

    l_i, e_true, lr_pos = make_soft_case(rng, mode, p=p)
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
            # e_hat are the BM/Chien roots (corrections relative to r_cand)
            e_hat_locs = roots[:]
            # total correction relative to original r_base:
            # e_total = flipmask XOR e_hat
            e_total_locs = xor_locations(flip_locs, e_hat_locs)

            # Build decoded codeword for correlation: c_hat = r_base XOR e_total
            c_hat = [rb for rb in r_base]
            for j in e_total_locs:
                c_hat[j] ^= 1
            corr = correlation(l_i, c_hat)
        else:
            e_hat_locs = []
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
            "roots": roots,              # e_hat locations (relative to r_cand)
            "e_total_locs": e_total_locs, # FINAL correction locations (relative to r_base)
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

    # answer: FINAL correction locations (flipmask XOR roots)
    with open(fans, "w") as f_ans:
        if len(win_total) == 0:
            f_ans.write(u10_bits(1023) + "\n")
        else:
            for loc in win_total:
                f_ans.write(u10_bits(loc) + "\n")

    # HR
    with open(fhr, "w") as f:
        f.write(f"Pattern p{pid} (SOFT, p=2, corrected answers include flipmask)\n")
        f.write(f"Mode code={mode.code}, (n,k)=({mode.n},{mode.k}), m={mode.m}, t={mode.t}\n")
        f.write(f"Primitive poly = 0x{mode.prim_poly:X}\n")
        f.write(f"Chien convention: root_at_inverse={root_at_inverse}, loc_is_degree={loc_is_degree}\n\n")

        f.write("1) Original channel LLRs l_i for r_i (i=0..n-1, ascending degree):\n")
        f.write("   i : l_i (hard_bit)\n")
        for i in range(mode.n):
            hb = 0 if l_i[i] >= 0 else 1
            f.write(f"   {i:4d}: {l_i[i]:4d} ({hb})\n")
        f.write("\n")
        f.write(f"True error locations (tx all-zero): {e_true}\n\n")

        f.write("2) Two least reliable indices:\n")
        f.write(f"   constructed lr_pos = {lr_pos}\n")
        f.write(f"   found by |LLR| sort = ({lrA},{lrB})\n")
        f.write(f"   lrA={lrA} l={l_i[lrA]} |l|={abs(l_i[lrA])}\n")
        f.write(f"   lrB={lrB} l={l_i[lrB]} |l|={abs(l_i[lrB])}\n\n")

        f.write("3/4/5) Chase variants (4 patterns):\n")
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="patterns_soft_p2_corrected", help="output directory")
    ap.add_argument("--num_patterns", type=int, default=3, help="patterns per mode (default: 3)")
    ap.add_argument("--seed", type=int, default=20251213, help="global seed")
    ap.add_argument("--chien_root_at_alpha_j", action="store_true", default=False,
                    help="Use sigma(alpha^j)==0 convention (default is sigma(alpha^{-j})==0)")
    ap.add_argument("--loc_is_msb", action="store_true", default=False,
                    help="Return location (n-1-j) instead of j")
    args = ap.parse_args()

    root_at_inverse = not args.chien_root_at_alpha_j
    loc_is_degree = not args.loc_is_msb

    base = {1: 400, 2: 500, 3: 600}
    for code, mode in MODES.items():
        for pidx in range(args.num_patterns):
            pid = base[code] + pidx
            patt_seed = (args.seed * 1000003) ^ (pid * 9176) ^ (code * 1315423911)
            write_one_soft_pattern(
                outdir=args.outdir,
                pid=pid,
                mode=mode,
                seed=patt_seed,
                p=2,
                root_at_inverse=root_at_inverse,
                loc_is_degree=loc_is_degree,
            )


if __name__ == "__main__":
    main()
