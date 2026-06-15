import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from base import set_state, out
from mod_add import mod_add, cmod_add
from add_const import add_const, cadd_const, ccadd_const
from doubling import doubling


def test_doubling():
    for i in range(7):
        c = tc.Circuit(5)
        c.append(set_state(i, [1, 2, 3, 4]))
        c.append(doubling([1, 2, 3, 4], [0]))
        assert out(c.state(), 1, 4) == (2 * i) % 7, f"test doubling with parameter {i} has problem, expected {(2 * i) % 7}, got {out(c.state(), 1, 4)}"
