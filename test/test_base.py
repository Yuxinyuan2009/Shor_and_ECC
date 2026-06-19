import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import tensorcircuit as tc
import numpy as np
from shor_ecc.base import set_state, out


def test_set_state():
    for i in range(16):
        c = tc.Circuit(4)
        c = set_state(i, [0, 1, 2, 3])
        assert out(c.state()) == i, f"setting {i} state has problem, expected {i} but got {out(c.state())}"
