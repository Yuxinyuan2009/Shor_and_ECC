import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from shor_ecc.base import set_state, out, ccphase
from shor_ecc.QFT import QFT
from shor_ecc.add import add, cadd
from shor_ecc.add_const import add_const, cadd_const, ccadd_const
from shor_ecc.mod_add import mod_add, cmod_add
from shor_ecc.negation import negation, cnegation


def test_negation():
    for i in range(7):
        c = tc.Circuit(5)
        c.append(set_state(i, [1, 2, 3, 4]))
        c.append(negation([1, 2, 3, 4], [0]))
        assert out(c.state(), 1, 4) == (-i) % 7, f"test negation with parameter {i} has problem, expected {(-i) % 7}, got {out(c.state(), 1, 4)}"


def test_cnegation():
    for p_val in [0, 1]:
        for i in range(7):
            c = tc.Circuit(6)
            c.append(set_state(i, [1, 2, 3, 4]))
            if p_val == 1:
                c.x(0)
            c.append(cnegation(0, [1, 2, 3, 4], [5]))
            expected = (-i) % 7 if p_val == 1 else i
            assert out(c.state(), 1, 4) == expected, f"test cnegation with control {p_val} and parameter {i} has problem, expected {expected}, got {out(c.state(), 1, 4)}"
