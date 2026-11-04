VERILATOR ?= verilator
NODE ?= node
.PHONY: lint test
lint:
	$(VERILATOR) --lint-only --timing -Wno-fatal --top-module bch rtl/bch.v
test:
	$(NODE) --test tests/*.test.js
