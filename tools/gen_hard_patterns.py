#!/usr/bin/env python3
"""
Hard-decision BCH pattern generator for:
  code=1: (63,51),   m=6,  t=2,  n=63,  input LLR count = 64  (LLR0 don't care)
  code=2: (255,239), m=8,  t=2,  n=255, input LLR count = 256
  code=3: (1023,983),m=10, t=4,  n=1023,input LLR count = 1024

Pattern file pXXX.txt:
  - Contains streaming idata lines.
  - Each line is 64 bits, representing 8 signed 8-bit LLRs (two's complement).
  - Bit layout in a line is: [63:56]=LLR0 ... [7:0]=LLR7.
  - For each codeword we output (n+1)/8 lines (because there are n+1 LLRs).
  - LLR ordering: LLR0 is don't-care (we force 0), LLR1 maps to X^(n-1), ...,
    LLR(n) maps to X^0.

Answer file pXXXa.txt:
  - One 10-bit unsigned integer per output cycle (binary string).
  - For each codeword:
      * if v=0 errors: output one line "1111111111" (1023)
      * else: output v lines for locations in ascending order

Human-readable pXXXhr.txt:
  - For each codeword: error locs, syndromes S1..S(2t), sigma(x) coefficients.

Notes:
  - We generate transmitted codeword = all-zero (valid because BCH is linear),
    so received bits are exactly the injected error pattern.
  - Primitive polynomials:
      m=6:  x^6 + x + 1
      m=8:  x^8 + x^4 + x^3 + x^2 + 1
      m=10: x^10 + x^3 + 1
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
    prim_poly: int  # includes x^m term (bit m set)
    exp: List[int]  # exp[k] = alpha^k
    log: List[int]  # log[val] = k, for val!=0

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
            red = prim_poly & mask  # lower m bits of modulus (excluding x^m)
            for _ in range(m * 2):  # safe upper bound
                if bb & 1:
                    res ^= aa
                bb >>= 1
                if bb == 0:
                    break
                carry = aa & (1 << (m - 1))
                aa = (aa << 1) & mask
                if carry:
                    aa ^= red
            return res & mask

        # alpha = x = 0b10 in polynomial basis
        exp[0] = 1
        for k in range(1, size - 1):
            exp[k] = mul(exp[k - 1], 2)
        for k, v in enumerate(exp):
            log[v] = k
        return GF2m(m=m, prim_poly=prim_poly, exp=exp, log=log)

    def add(self, a: int, b: int) -> int:
        return a ^ b

    def mul(self, a: int, b: int) -> int:
        # Russian peasant multiply with reduction by prim_poly
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

    def alpha_pow(self, k: int) -> int:
        if k < 0:
            k %= (len(self.exp))
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
        return self.n + 1  # includes LLR0 don't care

    @property
    def lines_per_codeword(self) -> int:
        assert self.llr_count % 8 == 0
        return self.llr_count // 8


MODES = {
    # (63,51) BCH: m=6, n=63, t=2, primitive poly x^6 + x + 1
    1: BCHMode(code=1, n=63, k=51, t=2, m=6, prim_poly=0x43),
    # (255,239) BCH: m=8, n=255, t=2, primitive poly x^8 + x^4 + x^3 + x^2 + 1
    2: BCHMode(code=2, n=255, k=239, t=2, m=8, prim_poly=0x11D),
    # (1023,983) BCH: m=10, n=1023, t=4, primitive poly x^10 + x^3 + 1
    3: BCHMode(code=3, n=1023, k=983, t=4, m=10, prim_poly=0x409),
}


# -------------------------
# Formatting helpers
# -------------------------
def int8_to_twos_comp_bits(x: int) -> str:
    """Convert signed int8 (-128..127) to 8-bit two's complement string."""
    if not (-128 <= x <= 127):
        raise ValueError(f"int8 out of range: {x}")
    if x < 0:
        x = (1 << 8) + x
    return format(x, "08b")


def pack_8_llrs_to_line(llrs: List[int]) -> str:
    """Pack 8 signed int8 LLRs into a 64-bit line string, LLR0 first (bits 63:56)."""
    assert len(llrs) == 8
    return "".join(int8_to_twos_comp_bits(v) for v in llrs)


def u10_bits(x: int) -> str:
    if not (0 <= x <= 1023):
        raise ValueError(f"u10 out of range: {x}")
    return format(x, "010b")


# -------------------------
# Core generation logic
# -------------------------
def choose_error_locations(rng: random.Random, n: int, t: int, allow_zero: bool = True) -> List[int]:
    v = rng.randint(0 if allow_zero else 1, t)
    if v == 0:
        return []
    locs = rng.sample(range(n), v)
    locs.sort()
    return locs


def build_received_bits_from_errors(n: int, error_locs: List[int]) -> List[int]:
    """
    Transmitted codeword is all-zero.
    Received r_j = 1 iff j in error_locs.
    r is indexed by polynomial degree j: r0 + r1 X + ... + r_{n-1} X^{n-1}.
    """
    r = [0] * n
    for j in error_locs:
        r[j] ^= 1
    return r


def build_llrs_for_stream(n: int, r_bits_deg_asc: List[int], rng: random.Random,
                          mag_min: int = 20, mag_max: int = 100) -> List[int]:
    """
    Returns LLR[0..n] where:
      - LLR0 is don't-care (forced 0)
      - LLR1 corresponds to r_{n-1} (X^{n-1})
      - ...
      - LLRn corresponds to r_0 (X^0)
    LLR sign rule: LLR >= 0 => bit tends to 0, LLR < 0 => bit tends to 1.
    """
    assert len(r_bits_deg_asc) == n
    llrs = [0] * (n + 1)
    llrs[0] = 0  # guaranteed 00000000 as you described
    for i in range(1, n + 1):
        bit = r_bits_deg_asc[n - i]  # i=1 -> r_{n-1}
        mag = rng.randint(mag_min, mag_max)
        llrs[i] = +mag if bit == 0 else -mag
    return llrs


def syndromes_from_error_locs(gf: GF2m, error_locs: List[int], t: int) -> List[int]:
    """Compute S1..S(2t) given e(X)=sum X^j and S_i = e(alpha^i) = sum alpha^(i*j)."""
    S = []
    for i in range(1, 2 * t + 1):
        acc = 0
        for j in error_locs:
            acc ^= gf.alpha_pow(i * j)
        S.append(acc)
    return S


def sigma_from_error_locs(gf: GF2m, error_locs: List[int]) -> List[int]:
    """
    Locator polynomial: sigma(x) = Π (1 + alpha^j * x)
    Coeff list is [sigma0, sigma1, ..., sigmad] for x^0..x^d.
    """
    sigma = [1]  # start with 1
    for j in error_locs:
        a_j = gf.alpha_pow(j)
        # multiply sigma by (1 + a_j x)
        new = [0] * (len(sigma) + 1)
        for d, c in enumerate(sigma):
            new[d] ^= c               # times 1
            new[d + 1] ^= gf.mul(c, a_j)  # times a_j x
        sigma = new
    return sigma


def write_pattern_files(outdir: str, pid: int, mode: BCHMode, num_codewords: int,
                        seed: int, mag_min: int, mag_max: int) -> None:
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)

    gf = GF2m.build(mode.m, mode.prim_poly)

    fname_in = os.path.join(outdir, f"p{pid}.txt")
    fname_ans = os.path.join(outdir, f"p{pid}a.txt")
    fname_hr = os.path.join(outdir, f"p{pid}hr.txt")

    with open(fname_in, "w") as f_in, open(fname_ans, "w") as f_ans, open(fname_hr, "w") as f_hr:
        f_hr.write(f"Pattern p{pid}\n")
        f_hr.write(f"Mode code={mode.code}, (n,k)=({mode.n},{mode.k}), m={mode.m}, t={mode.t}\n")
        f_hr.write(f"Primitive poly = 0x{mode.prim_poly:X}\n")
        f_hr.write(f"Codewords = {num_codewords}\n\n")

        for cw_idx in range(num_codewords):
            # choose <= t errors so hard-decision must succeed
            error_locs = choose_error_locations(rng, mode.n, mode.t, allow_zero=True)

            r_bits = build_received_bits_from_errors(mode.n, error_locs)
            llrs = build_llrs_for_stream(mode.n, r_bits, rng, mag_min=mag_min, mag_max=mag_max)

            # write input stream lines (8 LLRs per 64-bit line)
            for line_i in range(mode.lines_per_codeword):
                chunk = llrs[line_i * 8:(line_i + 1) * 8]
                f_in.write(pack_8_llrs_to_line(chunk) + "\n")

            # write answer stream (as odata would appear when finish is high)
            if len(error_locs) == 0:
                f_ans.write(u10_bits(1023) + "\n")
            else:
                for loc in error_locs:
                    f_ans.write(u10_bits(loc) + "\n")

            # human readable
            S = syndromes_from_error_locs(gf, error_locs, mode.t)
            sigma = sigma_from_error_locs(gf, error_locs)

            f_hr.write(f"=== Codeword {cw_idx:03d} ===\n")
            f_hr.write(f"errors v={len(error_locs)} locs={error_locs}\n")

            f_hr.write("Syndromes:\n")
            for i, s in enumerate(S, start=1):
                f_hr.write(
                    f"  S{i}: bits={gf.to_bits(s)}  alpha={gf.to_alpha_str(s)}  int={s}\n"
                )

            f_hr.write("Sigma(x)=Π(1 + a^loc * x):\n")
            for deg, c in enumerate(sigma):
                f_hr.write(
                    f"  sigma[{deg}] (x^{deg}): bits={gf.to_bits(c)}  alpha={gf.to_alpha_str(c)}  int={c}\n"
                )
            f_hr.write("\n")

    print(f"[OK] wrote {fname_in}, {fname_ans}, {fname_hr}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="patterns_hard", help="output directory")
    ap.add_argument("--num_codewords", type=int, default=100, help="codewords per pattern (default: 100)")
    ap.add_argument("--num_patterns", type=int, default=3, help="patterns per mode (default: 3)")
    ap.add_argument("--seed", type=int, default=12345, help="global seed (default: 12345)")
    ap.add_argument("--mag_min", type=int, default=20, help="min |LLR| for nonzero bits")
    ap.add_argument("--mag_max", type=int, default=100, help="max |LLR| for nonzero bits")
    args = ap.parse_args()

    # naming convention bases
    base = {1: 101, 2: 201, 3: 301}

    # generate patterns per mode
    for code, mode in MODES.items():
        for p in range(args.num_patterns):
            pid = base[code] + p
            # derive a per-pattern seed so outputs are repeatable
            patt_seed = (args.seed * 1000003) ^ (pid * 9176) ^ (code * 1315423911)
            write_pattern_files(
                outdir=args.outdir,
                pid=pid,
                mode=mode,
                num_codewords=args.num_codewords,
                seed=patt_seed,
                mag_min=args.mag_min,
                mag_max=args.mag_max,
            )


if __name__ == "__main__":
    main()
