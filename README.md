# Reconfigurable BCH Decoder

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
