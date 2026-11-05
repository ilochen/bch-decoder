#!/usr/bin/env python3
"""
Edge-case hard-decision BCH pattern generator.
Produces "difficult" codewords for each mode while ensuring #errors <= t.

Outputs (same format as previous generator):
  pXXX.txt   : packed 64-bit lines, each line contains 8 signed 8-bit LLRs
  pXXXa.txt  : 10-bit unsigned error locations (ascending) or 1023 if no errors
  pXXXhr.txt : human readable: error locs, syndromes, sigma coefficients

Pattern IDs used (customizable):
  m=6  -> p151, p152, ...
  m=8  -> p251, p252, ...
  m=10 -> p351, p352, ...
"""

from __future__ import annotations
import argparse
import os
import random
from dataclasses import dataclass
from typing import List, Tuple, Optional


# -------------------------
# GF(2^m) arithmetic
# -------------------------
@dataclass(frozen=True)
class GF2m:
    m: int
    prim_poly: int  # includes x^m term
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

    def alpha_pow(self, k: int) -> int:
        return self.exp[k % len(self.exp)]

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
    1: BCHMode(code=1, n=63,   k=51,  t=2, m=6,  prim_poly=0x43),   # x^6 + x + 1
    2: BCHMode(code=2, n=255,  k=239, t=2, m=8,  prim_poly=0x11D),  # x^8 + x^4 + x^3 + x^2 + 1
    3: BCHMode(code=3, n=1023, k=983, t=4, m=10, prim_poly=0x409),  # x^10 + x^3 + 1
}


# -------------------------
# Packing helpers
# -------------------------
def int8_to_twos_comp_bits(x: int) -> str:
    if not (-128 <= x <= 127):
        raise ValueError(f"int8 out of range: {x}")
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
# BCH math helpers
# -------------------------
def syndromes_from_error_locs(gf: GF2m, error_locs: List[int], t: int) -> List[int]:
    # S_i = sum alpha^(i*j) for each error location j
    S = []
    for i in range(1, 2 * t + 1):
        acc = 0
        for j in error_locs:
            acc ^= gf.alpha_pow(i * j)
        S.append(acc)
    return S


def sigma_from_error_locs(gf: GF2m, error_locs: List[int]) -> List[int]:
    # sigma(x) = Π (1 + alpha^j x)
    sigma = [1]
    for j in error_locs:
        a_j = gf.alpha_pow(j)
        new = [0] * (len(sigma) + 1)
        for d, c in enumerate(sigma):
            new[d] ^= c
            new[d + 1] ^= gf.mul(c, a_j)
        sigma = new
    return sigma


def build_received_bits(n: int, error_locs: List[int]) -> List[int]:
    r = [0] * n
    for j in error_locs:
        r[j] ^= 1
    return r


# -------------------------
# LLR synthesis (edge stressing)
# -------------------------
@dataclass(frozen=True)
class LLRStyle:
    name: str
    nonerr_mag: Tuple[int, int]   # range for non-error magnitudes
    err_mag: Tuple[int, int]      # range for error magnitudes (negative)
    zero_fraction_nonerr: float   # fraction of non-error bits set to 0 LLR
    saturate_fraction: float      # fraction of bits forced to +/-127/ -128

def make_llrs(mode: BCHMode, r_bits_deg_asc: List[int], rng: random.Random, style: LLRStyle) -> List[int]:
    """
    LLR[0..n], with mapping:
      LLR1 -> r_{n-1}, ... , LLRn -> r_0
    Rule: LLR >= 0 -> hard 0, LLR < 0 -> hard 1
    Keep LLR0 = 0.
    """
    n = mode.n
    llrs = [0] * (n + 1)
    llrs[0] = 0

    # Pre-decide saturation mask for positions 1..n
    sat_mask = [False] * (n + 1)
    if style.saturate_fraction > 0:
        for i in range(1, n + 1):
            if rng.random() < style.saturate_fraction:
                sat_mask[i] = True

    for i in range(1, n + 1):
        bit = r_bits_deg_asc[n - i]  # i=1 -> r_{n-1}
        if bit == 1:
            # error bit must be negative (cannot be 0)
            if sat_mask[i]:
                llrs[i] = -128
            else:
                mag = rng.randint(style.err_mag[0], style.err_mag[1])
                mag = max(1, mag)
                llrs[i] = -mag
        else:
            # non-error bit: positive, or sometimes 0 to create "tie"
            if sat_mask[i]:
                llrs[i] = 127
            else:
                if rng.random() < style.zero_fraction_nonerr:
                    llrs[i] = 0
                else:
                    mag = rng.randint(style.nonerr_mag[0], style.nonerr_mag[1])
                    llrs[i] = +mag

    return llrs


# -------------------------
# Edge-case builders
# -------------------------
def fixed_edge_error_sets(mode: BCHMode) -> List[List[int]]:
    n, t = mode.n, mode.t
    cases: List[List[int]] = []

    # Always include 0-error case
    cases.append([])

    # Single-error edges + midpoints
    cases += [[0], [n - 1], [1], [n - 2], [n // 2], [n // 2 - 1]]

    if t >= 2:
        # Two-error classic nasty cases
        cases += [
            [0, 1],
            [n - 2, n - 1],
            [0, n - 1],
            [1, n - 1],
            [0, n - 2],
            [n // 2 - 1, n // 2],
            [n // 2, n // 2 + 1],
            [n // 3, 2 * n // 3],
        ]

    if t >= 4:
        # Four-error nasty cases (bursts + split clusters + spread)
        cases += [
            [0, 1, 2, 3],
            [n - 4, n - 3, n - 2, n - 1],
            [0, 1, n - 2, n - 1],
            [0, 1, 2, n - 1],
            [0, n // 3, 2 * n // 3, n - 1],
            [n // 2 - 2, n // 2 - 1, n // 2, n // 2 + 1],
        ]

    # Ensure within [0,n-1] and unique/sorted
    norm = []
    for c in cases:
        cc = sorted(set([x for x in c if 0 <= x < n]))
        if len(cc) <= t:
            norm.append(cc)
    return norm


def find_case_with_zero_syndrome(gf: GF2m, n: int, t: int, target_Si: int,
                                max_tries: int, rng: random.Random) -> Optional[List[int]]:
    """
    Try to find error locations (<=t, >=1 errors) such that S_target_Si == 0.
    We bias toward using exactly t errors for harder structure.
    """
    for _ in range(max_tries):
        v = t if t > 1 else 1
        locs = rng.sample(range(n), v)
        locs.sort()
        S = syndromes_from_error_locs(gf, locs, t)
        if 1 <= target_Si <= len(S) and S[target_Si - 1] == 0:
            return locs
    return None


def find_case_with_zero_sigma_coeff(gf: GF2m, n: int, t: int,
                                   max_tries: int, rng: random.Random) -> Optional[List[int]]:
    """
    Only meaningful for t>=4: try to find <=t errors (use t errors) where sigma has a zero coefficient.
    """
    if t < 4:
        return None
    for _ in range(max_tries):
        locs = rng.sample(range(n), t)
        locs.sort()
        sigma = sigma_from_error_locs(gf, locs)
        # Ignore sigma[0]=1, look for any exact zero in middle/upper coeffs
        if any(c == 0 for c in sigma[1:]):
            return locs
    return None


def build_edge_case_suite(mode: BCHMode, seed: int, desired: int = 100) -> List[Tuple[List[int], LLRStyle]]:
    rng = random.Random(seed)
    gf = GF2m.build(mode.m, mode.prim_poly)

    # LLR styles that are "hard safe" but soft nasty
    styles = [
        LLRStyle("near_zero_errors", nonerr_mag=(40, 100), err_mag=(1, 3),  zero_fraction_nonerr=0.10, saturate_fraction=0.0),
        LLRStyle("many_zero_nonerr", nonerr_mag=(30, 80),  err_mag=(1, 5),  zero_fraction_nonerr=0.50, saturate_fraction=0.0),
        LLRStyle("saturation_mix",   nonerr_mag=(20, 60),  err_mag=(1, 10), zero_fraction_nonerr=0.10, saturate_fraction=0.15),
        LLRStyle("extreme_clean",    nonerr_mag=(90, 127), err_mag=(90, 128),zero_fraction_nonerr=0.0,  saturate_fraction=0.25),
    ]

    cases: List[Tuple[List[int], LLRStyle]] = []
    base_sets = fixed_edge_error_sets(mode)

    # Pair each fixed error set with a rotating style
    for i, locs in enumerate(base_sets):
        cases.append((locs, styles[i % len(styles)]))

    # Add automatically found "structural" cases
    # Try to find some with zero syndrome at various indices
    for target_Si in range(2, 2 * mode.t + 1):
        locs = find_case_with_zero_syndrome(gf, mode.n, mode.t, target_Si=target_Si,
                                            max_tries=20000 if mode.n <= 255 else 60000, rng=rng)
        if locs is not None:
            cases.append((locs, styles[target_Si % len(styles)]))

    # For t=4 mode, try to find sigma with zero coefficients
    locs_sigma0 = find_case_with_zero_sigma_coeff(gf, mode.n, mode.t,
                                                  max_tries=120000, rng=rng)
    if locs_sigma0 is not None:
        cases.append((locs_sigma0, styles[1]))

    # Fill up to desired count by generating rotated/consecutive/spread variants
    while len(cases) < desired:
        pick_locs, pick_style = cases[rng.randrange(len(cases))]
        v = len(pick_locs)

        if v == 0:
            # keep some all-zero / trivial cases but vary LLR style
            cases.append(([], styles[rng.randrange(len(styles))]))
            continue

        # Create a variant by shifting all locations by a random offset mod n
        shift = rng.randrange(mode.n)
        new_locs = sorted({(x + shift) % mode.n for x in pick_locs})
        # If collisions reduced v too much, fix by adding extra random locs up to v
        while len(new_locs) < v:
            new_locs.append(rng.randrange(mode.n))
            new_locs = sorted(set(new_locs))

        # Keep within t (it is), and prefer edge-ish: sometimes force include 0 or n-1
        if rng.random() < 0.25 and 0 not in new_locs:
            new_locs[0] = 0
            new_locs = sorted(set(new_locs))
        if rng.random() < 0.25 and (mode.n - 1) not in new_locs:
            new_locs[-1] = mode.n - 1
            new_locs = sorted(set(new_locs))

        # If we accidentally exceeded t due to set operations, trim deterministically
        if len(new_locs) > mode.t:
            new_locs = new_locs[:mode.t]

        cases.append((new_locs, styles[rng.randrange(len(styles))]))

    return cases[:desired]


# -------------------------
# File writer
# -------------------------
def write_edge_pattern(outdir: str, pid: int, mode: BCHMode, suite: List[Tuple[List[int], LLRStyle]], seed: int) -> None:
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    gf = GF2m.build(mode.m, mode.prim_poly)

    fin = os.path.join(outdir, f"p{pid}.txt")
    fans = os.path.join(outdir, f"p{pid}a.txt")
    fhr = os.path.join(outdir, f"p{pid}hr.txt")

    with open(fin, "w") as f_in, open(fans, "w") as f_ans, open(fhr, "w") as f_hr:
        f_hr.write(f"Pattern p{pid} (EDGE CASES)\n")
        f_hr.write(f"Mode code={mode.code}, (n,k)=({mode.n},{mode.k}), m={mode.m}, t={mode.t}\n")
        f_hr.write(f"Primitive poly = 0x{mode.prim_poly:X}\n")
        f_hr.write(f"Codewords = {len(suite)}\n\n")

        for cw_idx, (error_locs, style) in enumerate(suite):
            r_bits = build_received_bits(mode.n, error_locs)
            llrs = make_llrs(mode, r_bits, rng, style)

            # write packed idata lines
            for li in range(mode.lines_per_codeword):
                chunk = llrs[li * 8:(li + 1) * 8]
                f_in.write(pack_8_llrs_to_line(chunk) + "\n")

            # answer file: locations ascending, or 1023 once if no errors
            if len(error_locs) == 0:
                f_ans.write(u10_bits(1023) + "\n")
            else:
                for loc in error_locs:
                    f_ans.write(u10_bits(loc) + "\n")

            # human readable
            S = syndromes_from_error_locs(gf, error_locs, mode.t)
            sigma = sigma_from_error_locs(gf, error_locs)

            f_hr.write(f"=== Codeword {cw_idx:03d} ===\n")
            f_hr.write(f"LLR_style={style.name}\n")
            f_hr.write(f"errors v={len(error_locs)} locs={error_locs}\n")

            f_hr.write("Syndromes:\n")
            for i, s in enumerate(S, start=1):
                f_hr.write(f"  S{i}: bits={gf.to_bits(s)}  alpha={gf.to_alpha_str(s)}  int={s}\n")

            f_hr.write("Sigma(x)=Π(1 + a^loc * x):\n")
            for deg, c in enumerate(sigma):
                f_hr.write(f"  sigma[{deg}] (x^{deg}): bits={gf.to_bits(c)}  alpha={gf.to_alpha_str(c)}  int={c}\n")
            f_hr.write("\n")

    print(f"[OK] wrote {fin}, {fans}, {fhr}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="patterns_hard_edge", help="output directory")
    ap.add_argument("--num_patterns", type=int, default=1, help="patterns per mode (default: 1)")
    ap.add_argument("--num_codewords", type=int, default=100, help="codewords per pattern (default: 100)")
    ap.add_argument("--seed", type=int, default=777, help="global seed (default: 777)")
    args = ap.parse_args()

    # Start IDs for edge patterns (so they won't collide with your basic patterns)
    base = {1: 151, 2: 251, 3: 351}

    for code, mode in MODES.items():
        for p in range(args.num_patterns):
            pid = base[code] + p
            suite_seed = (args.seed * 1000003) ^ (pid * 9176) ^ (code * 1315423911)
            suite = build_edge_case_suite(mode, seed=suite_seed, desired=args.num_codewords)

            write_seed = suite_seed ^ 0xA5A5A5A5
            write_edge_pattern(args.outdir, pid, mode, suite, seed=write_seed)


if __name__ == "__main__":
    main()
