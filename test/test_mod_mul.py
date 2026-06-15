import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from base import set_state, out, ccphase
from QFT import QFT
from add import add, cadd
from add_const import add_const, cadd_const, ccadd_const
from mod_add import mod_add, cmod_add
from doubling import doubling
from mod_add_const import mod_add_const, cmod_add_const, ccmod_add_const, cccmod_add_const
from mod_mul import mod_mul, mod_square, cmod_square, inv


def test_mod_mul():
    for i in np.arange(7):
        for j in range(7):
            c = tc.Circuit(12)
            c.append(set_state(i, [0, 1, 2, 3]))
            c.append(set_state(j, [4, 5, 6, 7]))
            c.append(mod_mul([0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]))
            assert out(c.state(), 8, 11) == (i * j) % 7, f"test mod_mul with parameters {i}, {j} has problem, expected {(i * j) % 7}, got {out(c.state(), 8, 11)}"


def test_mod_square():
    for i in np.arange(7):
        for j in range(7):
            c = tc.Circuit(8)
            c.append(set_state(i, [0, 1, 2, 3]))
            c.append(set_state(j, [4, 5, 6, 7]))
            c.append(mod_square([0, 1, 2, 3], [4, 5, 6, 7]))
            assert out(c.state(), 4, 7) == (j + i * i) % 7, f"test mod_square with parameter {i}, {j} has problem, expected {(i * i) % 7}, got {out(c.state(), 4, 7)}"


def test_cmod_square():
    for p_val in [0, 1]:
        for i in np.arange(7):
            j = 3
            c = tc.Circuit(9)
            c.append(set_state(i, [1, 2, 3, 4]))
            c.append(set_state(j, [5, 6, 7, 8]))
            if p_val == 1:
                c.x(0)
            c.append(cmod_square(0, [1, 2, 3, 4], [5, 6, 7, 8]))
            expected = (j + i * i) % 7 if p_val == 1 else j
            assert out(c.state(), 5, 8) == expected, f"test cmod_square with parameters {p_val}, {i}, {j} has problem, expected {expected}, got {out(c.state(), 5, 8)}"


def test_inv():
    for i in range(7):
        c = tc.Circuit(4)
        c.append(set_state(i, [0, 1, 2, 3]))
        c.append(inv([0, 1, 2, 3]))
        assert out(c.state(), 0, 3) == (i**5) % 7, f"test inv with parameter {i} has problem, expected {(i**5) % 7}, got {out(c.state(), 0, 3)}"
