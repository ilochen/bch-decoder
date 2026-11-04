# Validation

The original RTL is preserved as supplied. `make lint` elaborates the design with Verilator; `-Wno-fatal` retains compiler warnings in the output.

Observed warning classes: CASEINCOMPLETE, MULTIDRIVEN, UNOPTFLAT, WIDTHEXPAND, WIDTHTRUNC.

`npm test` checks the accompanying architectural reference models against independent numerical or sequential references. These tests validate model behavior; RTL simulation and synthesis are separate flows.
