import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from base import set_state, out
from QFT import QFT
from add import add, cadd


def test_add():
    for i in range(16):
        for j in range(16):
            c = tc.Circuit(8)
            c.append(set_state(i, [0, 1, 2, 3]))
            c.append(set_state(j, [4, 5, 6, 7]))
            c.append(add([0, 1, 2, 3], [4, 5, 6, 7]))
            assert out(c.state(), 4, 7) == (i + j) % 16, f"test add with parameter {i}, {j} has problem, expected {(i + j) % 16}, got {out(c.state(), 4, 7)}"


def test_cadd():
    for p_val in [0, 1]:
        for i in range(16):
            for j in range(16):
                c = tc.Circuit(9)
                if p_val == 1:
                    c.X(8)
                c.append(set_state(i, [0, 1, 2, 3]))
                c.append(set_state(j, [4, 5, 6, 7]))
                c.append(cadd(8, [0, 1, 2, 3], [4, 5, 6, 7]))
                expected = (j + p_val * i) % 16
                assert out(c.state(), 4, 7) == expected, f"test cadd with parameters {p_val}, {i}, {j} has problem, expected {expected}, got {out(c.state(), 4, 7)}"
