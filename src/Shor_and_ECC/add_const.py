import tensorcircuit as tc
import numpy as np
from Shor_and_ECC.add import *       # Provides QFT and ccphase

@ block
def add_const(x: list, t: int) -> tc.Circuit:
    """
    Add a classical constant to a quantum register.

    |x> -> |x + c (mod 2**n)>
    where c is the constant integer `t` and n = len(x).

    Parameters
    ----------
    x : list
        Qubit indices for the n-qubit register (holds |x>).
    t : int
        Classical constant to add (must be < 2**n).

    Returns
    -------
    tc.Circuit
        Circuit implementing |x> → |x + t (mod 2^n)>.
    """
    c = tc.Circuit(max(x) + 1)
    n = len(x)

    # Almost copy the process in add, but replace the control qubit with the classical constant t
    c.append(QFT(x))
    for j in range(n):
        if (t >> (n - 1 - j)) & 1:
            for i in range(n - 1 - j, n):
                c.phase(x[i], theta=np.pi / (2 ** (i + j - n + 1)))
    c.append(QFT(x).inverse())
    return c


@ block
def cadd_const(p: int, x: list, t: int) -> tc.Circuit:
    """
    Controlled addition of a classical constant.

    |p>|x> -> |p>|x + p*c (mod 2**n)>
    where c = t (the constant). If control qubit p is |0>, nothing changes;
    if p is |1>, the constant is added.

    Parameters
    ----------
    p : int
        Index of the control qubit.
    x : list
        Qubit indices for the n-qubit target register.
    t : int
        Classical constant to add (when control is active).

    Returns
    -------
    tc.Circuit
        Circuit implementing the controlled constant addition.
    """
    c = tc.Circuit(max(x + [p]) + 1)
    n = len(x)
    c.append(QFT(x))
    for j in range(n):
        if (t >> (n - 1 - j)) & 1:
            for i in range(n - 1 - j, n):
                c.cphase(p, x[i], theta=np.pi / (2 ** (i + j - n + 1)))

    # Step 3: Inverse QFT
    c.append(QFT(x).inverse())
    return c


@ block
def ccadd_const(p: int, q: int, x: list, t: int) -> tc.Circuit:
    """
    Doubly-controlled addition of a classical constant.

    |p>|q>|x> -> |p>|q>|x + p*q*c (mod 2**n)>
    where c = t. The constant is added only when both control qubits p and q are |1>.

    Parameters
    ----------
    p : int
        Index of the first control qubit.
    q : int
        Index of the second control qubit.
    x : list
        Qubit indices for the n-qubit target register.
    t : int
        Classical constant to add (when both controls are active).

    Returns
    -------
    tc.Circuit
        Circuit implementing the doubly-controlled constant addition.
    """
    c = tc.Circuit(max(x + [p, q]) + 1)
    n = len(x)
    c.append(QFT(x))
    for j in range(n):
        if (t >> (n - 1 - j)) & 1:
            for i in range(n - 1 - j, n):
                c.append(ccphase(p, q, x[i], theta=np.pi / (2 ** (i + j - n + 1))))

    # Step 3: Inverse QFT
    c.append(QFT(x).inverse())
    return c

