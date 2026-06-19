import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from shor_ecc.base import set_state, out, ccphase
from shor_ecc.QFT import QFT
from shor_ecc.add import add
from shor_ecc.add_const import add_const, cadd_const, ccadd_const
from shor_ecc.mod_add_const import mod_add_const, cmod_add_const, ccmod_add_const


def test_mod_add_const():
    for i in range(7):
        for j in range(7):
            c = tc.Circuit(5)
            c.append(set_state(i, [1, 2, 3, 4]))
            c.append(mod_add_const([1, 2, 3, 4], [0], j))
            assert out(c.state(), 1, 4) == (i + j) % 7, f"test mod_add_const with parameters {i}, {j} has problem, expected {(i + j) % 7}, got {out(c.state(), 1, 4)}"


def test_cmod_add_const():
    for p_val in [0, 1]:
        for i in range(7):
            for j in range(7):
                c = tc.Circuit(6)
                c.append(set_state(i, [1, 2, 3, 4]))
                if p_val == 1:
                    c.x(0)
                c.append(cmod_add_const(0, [1, 2, 3, 4], [5], j))
                expected = (i + j) % 7 if p_val == 1 else i
                assert out(c.state(), 1, 4) == expected, f"test cmod_add_const with parameters {p_val}, {i}, {j} has problem, expected {expected}, got {out(c.state(), 1, 4)}"


def test_ccmod_add_const():
    for p_val in [0, 1]:
        for q_val in [0, 1]:
            for i in range(7):
                for j in range(7):
                    c = tc.Circuit(7)
                    c.append(set_state(i, [2, 3, 4, 5]))
                    if p_val == 1:
                        c.x(0)
                    if q_val == 1:
                        c.x(1)
                    c.append(ccmod_add_const(0, 1, [2, 3, 4, 5], [6], j))
                    expected = (i + j) % 7 if p_val == 1 and q_val == 1 else i
                    assert out(c.state(), 2, 5) == expected, f"test ccmod_add_const with parameters {p_val}, {q_val}, {i}, {j} has problem, expected {expected}, got {out(c.state(), 2, 5)}"
