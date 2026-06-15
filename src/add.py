import tensorcircuit as tc
import numpy as np
from QFT import *          # Provides QFT (Quantum Fourier Transform)
from base import *        # Provides set_state, out, ccphase (base utilities)

@ block
def add(x: list, y: list) -> tc.Circuit:
    """
    Quantum adder: |x>|y> -> |x>|x+y (mod 2**n)>
    Both x and y are lists of qubit indices representing n-qubit registers.

    Parameters
    ----------
    x : list
        Qubit indices for the first register (holds |x>).
    y : list
        Qubit indices for the second register (holds |y>, becomes |x+y>).

    Returns
    -------
    tc.Circuit
        Circuit implementing the in-place addition.
    """
    c = tc.Circuit(max(x + y) + 1)
    n = len(y)
    c.append(QFT(y))
    for i in range(n):
        for j in range(n - 1 - i, n):
            c.cphase(x[i], y[j], theta=np.pi / (2 ** (i + j - n + 1)))
    c.append(QFT(y).inverse())
    return c


@ block
def cadd(p: int, x: list, y: list) -> tc.Circuit:
    """
    Controlled quantum adder: |p>|x>|y> -> |p>|x>|y + p*x (mod 2**n)>
    The single control qubit p decides whether the addition is performed.
    If p = |1>, the adder acts normally; if p = |0>, nothing changes.

    Parameters
    ----------
    p : int
        Index of the control qubit.
    x : list
        Qubit indices of the first register (|x>).
    y : list
        Qubit indices of the second register (|y>, becomes |y + p*x>).

    Returns
    -------
    tc.Circuit
        Circuit implementing the controlled in-place addition.
    """
    c = tc.Circuit(max([p] + x + y) + 1)
    n = len(x)
    c.append(QFT(y))
    for i in range(n):
        for j in range(n - 1 - i, n):
            c.append(ccphase(p, x[i], y[j], theta=np.pi / (2 ** (i + j - n + 1))))
    c.append(QFT(y).inverse())
    return c


