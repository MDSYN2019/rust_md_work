# Simple runner targets for serial and MPI examples

CARGO ?= cargo
MPIRUN ?= mpirun
NP ?= 4

.PHONY: help run run-mpi build build-mpi check fmt test clean

help:
	@echo "Targets:"
	@echo "  make run        - Run default (non-MPI) MD example"
	@echo "  make run-mpi    - Run MPI-enabled MD example (set NP=<ranks>)"
	@echo "  make build      - Build default binary"
	@echo "  make build-mpi  - Build binary with MPI feature"
	@echo "  make check      - cargo check"
	@echo "  make fmt        - cargo fmt"
	@echo "  make test       - cargo test"
	@echo "  make clean      - cargo clean"

run:
	$(CARGO) run

run-mpi:
	$(MPIRUN) -n $(NP) $(CARGO) run --features mpi

build:
	$(CARGO) build

build-mpi:
	$(CARGO) build --features mpi

check:
	$(CARGO) check

fmt:
	$(CARGO) fmt

test:
	$(CARGO) test

clean:
	$(CARGO) clean
