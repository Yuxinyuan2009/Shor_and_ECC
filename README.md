# Shor's Algorithm for Elliptic Curve Cryptography

A quantum circuit simulation of **Shor's algorithm applied to the Elliptic Curve Discrete Logarithm Problem (ECDLP)**. Uses the `tensorcircuit` library to construct and simulate the quantum gates needed for elliptic curve point addition over GF(7).

## Overview

This project demonstrates how Shor's algorithm can break ECC by solving the ECDLP on a quantum computer simulator. It implements the full stack of quantum arithmetic operations needed for elliptic curve point addition:

- **Finite field arithmetic** over GF(7): addition, multiplication, squaring, negation, doubling
- **Elliptic curve point addition** handling all corner cases (identity, negation, doubling, -2P)
- **Quantum Fourier Transform** (QFT) and its inverse
- **Shor's algorithm circuit** with controlled point additions and measurement
- **Post-processing** for key recovery using continued fractions

## Project Structure

```
src/
└── shor_ecc/
    ├── __init__.py              # Package init
    ├── base.py                  # Core utilities: set_state, out, ccphase, block decorator
    ├── QFT.py                   # Quantum Fourier Transform
    ├── add.py                   # Quantum adder: |x>|y> -> |x>|x+y>
    ├── add_const.py             # Addition of classical constant
    ├── mod_add.py               # Modular addition (mod 7)
    ├── mod_add_const.py         # Modular addition of constant (mod 7)
    ├── mod_mul.py               # Modular multiplication and squaring (mod 7)
    ├── doubling.py              # Doubling: |x> -> |2x mod 7>
    ├── negation.py              # Modular negation: |x> -> |-x mod 7>
    ├── ECC_base.py              # Elliptic curve point addition circuits
    ├── ECC_circuit.py           # Main Shor's algorithm circuit (22 qubits, ~110k gates)
    └── ECC_post_measurement.py  # Post-processing and private key recovery

test/                        # pytest test suite for all modules
assets/                      # Simulation output and visualizations
```

## Requirements

- Python == 3.12
- tensorcircuit >= 0.11.0
- numpy >= 1.26.4
- sympy >= 1.14.0
- tqdm >= 4.67.1 (development)

## Installation

```bash
pip install .
```

For development:

```bash
pip install -e ".[dev]"
```

## Usage

Run the full Shor's algorithm simulation:

```bash
export PYTHONPATH="$PWD/src"
python -m shor_ecc.ECC_circuit
```

This outputs probability amplitudes to `assets/output.json` and prints runtime (approx. 197 seconds for 10 control qubits).

Run post-processing to recover the private key:

```bash
export PYTHONPATH="$PWD/src"
python -m shor_ecc.ECC_post_measurement
```

Run tests:

```bash
pytest
```

## How It Works

The circuit uses 22 qubits with the curve `y^2 = x^3 + ax + b (mod 7)`:

1. **Superposition**: Hadamard gates create superposition over 10 control qubits
2. **Controlled point additions**: Iteratively apply `cond_ECC_add`, doubling the point each step
3. **Inverse QFT**: Transform the period information into phase measurements
4. **Key recovery**: Post-process measurement results to extract the private key

## Results

- Circuit: 22 qubits, ~110,000 gates
- Runtime: ~197 seconds
- Elliptic curve: GF(7), `y^2 = x^3 + ax + b (mod 7)`