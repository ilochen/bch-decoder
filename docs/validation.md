# Validation

The original RTL is preserved as supplied. `make lint` elaborates the design with Verilator; `-Wno-fatal` retains compiler warnings in the output.

Observed warning classes: CASEINCOMPLETE, MULTIDRIVEN, UNOPTFLAT, WIDTHEXPAND, WIDTHTRUNC.

`npm test` checks the accompanying architectural reference models against independent numerical or sequential references. These tests validate model behavior; RTL simulation and synthesis are separate flows.

The Verilator RTL testbench passed fixtures 100, 200, 300, 400, 500 and 600 on Verilator 5.052, checking 27 output positions in total. Input and expected-answer files were selected from the supplied final-project test vectors. No cell library is required for this functional run.
