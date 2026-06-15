import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from base import set_state, out, ccphase
from QFT import QFT
from add import add, cadd
from add_const import add_const, cadd_const, ccadd_const
from mod_add import mod_add, cmod_add


def test_mod_add():
    for i in range(7):
        for j in range(7):
            c = tc.Circuit(9)
            c.append(set_state(i, [0, 1, 2, 3]))
            c.append(set_state(j, [4, 5, 6, 7]))
            c.append(mod_add(list(range(4)), list(range(4, 8)), [8]))
            assert out(c.state(), 4, 7) == (i + j) % 7, f"test mod_add with parameters {i}, {j} has problem"


def test_cmod_add():
    for p_val in [0, 1]:
        for i in range(7):
            for j in range(7):
                c = tc.Circuit(10)
                c.append(set_state(i, [1, 2, 3, 4]))
                c.append(set_state(j, [5, 6, 7, 8]))
                if p_val == 1:
                    c.x(0)
                c.append(cmod_add(0, list(range(1, 5)), list(range(5, 9)), [9]))
                expected = (i + j) % 7 if p_val == 1 else j
                assert out(c.state(), 5, 8) == expected, f"test cmod_add with parameters {p_val}, {i}, {j} has problem"
