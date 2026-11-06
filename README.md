# Reconfigurable BCH Decoder

[![Validation](https://github.com/ilochen/bch-decoder/actions/workflows/check.yml/badge.svg)](https://github.com/ilochen/bch-decoder/actions/workflows/check.yml)

A configurable hard/soft BCH decoder for BCH(63,51), BCH(255,239), and BCH(1023,983). The datapath combines byte-wide syndrome accumulation, inversionless Berlekamp–Massey, and eight-way Chien search. Four Chase candidates reuse syndrome computations; a two-phase schedule shares GF multipliers, and mode-dependent clock gating controls unused candidate state.

[Interactive architecture and design decisions](https://ilochen.github.io/projects/bch-decoder/)

## Source layout

- `rtl/bch.v`: complete decoder RTL, top module `bch`.
- `assets/bch-hardware.js`: architectural trace model for streaming syndrome, scaled BM, and eight-lane Chien search.
- `assets/models.js`: mathematical reference routines, including finite-field and BCH operations.
- `tests/bch-hardware.test.js`: independent polynomial-evaluation, root-location, bypass, and soft-candidate checks across all three code lengths.

## Run

```sh
make lint          # Verilator source elaboration
npm test           # Node.js 22+, no dependencies
```

The RTL interface includes `clk`, active-low `rstn`, mode/code selection, `set`, 64-bit `idata`, and `ready`/`finish` handshakes with 10-bit error-position output. Field elements use polynomial-basis bit representations. Constant multiplications reduce to XOR maps; variable GF multiplication uses AND partial products and modular reduction.

See [validation notes](docs/validation.md) for the checked artifacts and compiler diagnostics.

## Decode a fixture through the RTL

Run `make rtl-test` with Verilator 5 and a C++20 compiler. This builds the submitted `bch` module and drives its `set/ready/finish` interface with binary fixtures for **all three lengths, in both modes**. The checker compares each reported location, the number of responses, and the number of consumed input words. A mismatch or timeout exits with an error.

| Fixture | Code | Mode | Codewords | Input words | Output positions |
|---|---|---|---:|---:|---:|
| 100 | (63,51) | hard | 2 | 16 | 4 |
| 200 | (255,239) | hard | 2 | 64 | 4 |
| 300 | (1023,983) | hard | 2 | 256 | 8 |
| 400 | (63,51) | Chase | 1 | 8 | 2 |
| 500 | (255,239) | Chase | 1 | 32 | 4 |
| 600 | (1023,983) | Chase | 1 | 128 | 5 |

Each line in `sim/vectors/pNNN.txt` holds eight packed channel values (64 bits). Its `pNNNa.txt` partner holds ten-bit expected error positions in output order. Use `./build/Vdecoder_vectors +PATTERN=400` to replay one test after building. The six fixtures are smoke tests; the JavaScript tests independently check the polynomial and architectural models.

The four original Python generators in `tools/` generate hard/soft inputs and edge cases using finite-field reference arithmetic. Inspect each script's `--help` for its seed, case count and output directory. Generate experiments outside `sim/vectors/` to preserve the shipped regression inputs, for example:

```sh
python3 tools/gen_soft_patterns.py --num_patterns 1 --outdir build/generated-soft
```

## Implementation archive

- [Synthesis](implementation/synthesis/) contains the submitted mapped Verilog and the design-specific clock and interface constraints. The clock target is 5 ns. Configure libraries and operating conditions in the authorized synthesis environment.
- [Post-route](implementation/postroute/) contains the connectivity netlist exported after routing. The mapped and post-route netlists reference external standard-cell modules. Select one netlist **instead of** `rtl/bch.v` when simulating in an environment with authorized cell models.
- The submitted synthesis report gives a cell area of **587,717.96** in the tool's area units. The archived post-route timing summaries report **+0.563 ns** setup slack and **+1.245 ns** hold slack for the all-path groups, with zero violating paths in both summaries. These are tool estimates from the submitted run.

The [timing archive](docs/timing-archive.md) identifies the paired design-specific SDF exports, their hashes, and the requirements for a separate gate-level run.

## Validation snapshot

All six RTL fixtures passed: three code lengths in hard and Chase modes. Independent finite-field and architectural tests also pass.

Rechecked on 6 November 2026. Simulation results establish the covered functional behavior; physical implementation measurements are identified separately in the documentation.
