import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from shor_ecc.base import set_state, out, ccphase
from shor_ecc.QFT import QFT
from shor_ecc.add import add
from shor_ecc.add_const import add_const, cadd_const, ccadd_const


def test_add_const():
    for i in range(16):
        c = tc.Circuit(4)
        c.append(set_state(i, [0, 1, 2, 3]))
        c.append(add_const([0, 1, 2, 3], 5))
        assert out(c.state(), 0, 3) == (i + 5) % 16, f"test add_const with parameter {i} has problem, expected {(i + 5) % 16}, got {out(c.state(), 0, 3)}"


def test_cadd_const():
    for i in range(16):
        c = tc.Circuit(5)
        c.append(set_state(i, [1, 2, 3, 4]))
        c.x(0)
        c.append(cadd_const(0, [1, 2, 3, 4], 2))
        assert out(c.state(), 1, 4) == (i + 2) % 16, f"test cadd_const with parameter {i} has problem, expected {(i + 2) % 16}, got {out(c.state(), 1, 4)}"

    for i in range(16):
        c = tc.Circuit(5)
        c.append(set_state(i, [1, 2, 3, 4]))
        c.append(cadd_const(0, [1, 2, 3, 4], 2))
        assert out(c.state(), 1, 4) == i, f"test cadd_const with parameter {i} has problem, expected {i}, got {out(c.state(), 1, 4)}"


def test_ccadd_const():
    for i in range(16):
        c = tc.Circuit(6)
        c.append(set_state(i, [2, 3, 4, 5]))
        c.x(0)
        c.x(1)
        c.append(ccadd_const(0, 1, [2, 3, 4, 5], 3))
        assert out(c.state(), 2, 5) == (i + 3) % 16, f"test ccadd_const with parameter {i} has problem, expected {(i + 3) % 16}, got {out(c.state(), 2, 5)}"

    for i in range(16):
        c = tc.Circuit(6)
        c.append(set_state(i, [2, 3, 4, 5]))
        c.x(0)
        c.append(ccadd_const(0, 1, [2, 3, 4, 5], 3))
        assert out(c.state(), 2, 5) == i, f"test ccadd_const with parameter {i} has problem, expected {i}, got {out(c.state(), 2, 5)}"
