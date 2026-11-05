VERILATOR ?= verilator
NODE ?= node
.PHONY: lint test
lint:
	$(VERILATOR) --lint-only --timing -Wno-fatal --top-module bch rtl/bch.v
test:
	$(NODE) --test tests/*.test.js

.PHONY: rtl-test
rtl-test:
	mkdir -p build
	$(VERILATOR) --binary --timing -j 4 -Wno-fatal --Mdir build --top-module decoder_vectors sim/decoder_vectors.sv rtl/bch.v
	set -e; for p in 100 200 300 400 500 600; do ./build/Vdecoder_vectors +PATTERN=$$p; done
