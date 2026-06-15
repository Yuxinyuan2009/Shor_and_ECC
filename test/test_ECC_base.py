import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from base import set_state, out, cset_state, ccphase
from QFT import QFT
from add import add, cadd
from add_const import add_const, cadd_const, ccadd_const
from mod_add import mod_add, cmod_add
from doubling import doubling
from mod_add_const import mod_add_const, cmod_add_const, ccmod_add_const, cccmod_add_const
from mod_mul import mod_mul, mod_square, cmod_square, inv
from negation import negation, cnegation
from ECC_base import classical_ECC_add, ECC_add, cond_ECC_add


def test_ECC_add():
    c = tc.Circuit(12)
    c.append(set_state(3, [0, 1, 2, 3]))
    c.append(set_state(3, [4, 5, 6, 7]))
    c.append(ECC_add(5, 6, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == classical_ECC_add(3, 3, 5, 6, 2)


def test_cond_ECC_add_Q_eq_P():
    c = tc.Circuit(13)
    c.append(set_state(3, [0, 1, 2, 3]))
    c.append(set_state(3, [4, 5, 6, 7]))
    c.x(12)
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == classical_ECC_add(3, 3, 3, 3, 2)

    c = tc.Circuit(13)
    c.append(set_state(3, [0, 1, 2, 3]))
    c.append(set_state(3, [4, 5, 6, 7]))
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == (3, 3)


def test_cond_ECC_add_Q_eq_O():
    c = tc.Circuit(13)
    c.x(12)
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == classical_ECC_add(3, 3, 0, 0, 2)

    c = tc.Circuit(13)
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == (0, 0)


def test_cond_ECC_add_Q_eq_neg_P():
    c = tc.Circuit(13)
    c.append(set_state(3, [0, 1, 2, 3]))
    c.append(set_state(4, [4, 5, 6, 7]))
    c.x(12)
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == classical_ECC_add(3, 3, 3, 4, 2)

    c = tc.Circuit(13)
    c.append(set_state(3, [0, 1, 2, 3]))
    c.append(set_state(4, [4, 5, 6, 7]))
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == (3, 4)


def test_cond_ECC_add_Q_eq_neg_2P():
    c = tc.Circuit(13)
    c.append(set_state(2, [0, 1, 2, 3]))
    c.append(set_state(4, [4, 5, 6, 7]))
    c.x(12)
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == classical_ECC_add(3, 3, 2, 4, 2), f"{out(c.state(), 0, 3)}, {out(c.state(), 4, 7)}"

    c = tc.Circuit(13)
    c.append(set_state(2, [0, 1, 2, 3]))
    c.append(set_state(4, [4, 5, 6, 7]))
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == (2, 4)


def test_cond_ECC_add_normal():
    c = tc.Circuit(13)
    c.append(set_state(1, [0, 1, 2, 3]))
    c.append(set_state(0, [4, 5, 6, 7]))
    c.x(12)
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == classical_ECC_add(3, 3, 1, 0, 2)

    c = tc.Circuit(13)
    c.append(set_state(1, [0, 1, 2, 3]))
    c.append(set_state(0, [4, 5, 6, 7]))
    c.append(cond_ECC_add(12, 3, 3, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
    assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == (1, 0)
