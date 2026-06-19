import tensorcircuit as tc
import numpy as np
from typing import Union

def block(func):
	def wrapper(*args, **kwargs):
		c = func(*args, **kwargs)
		n = c._nqubits
		a = [i for i in range(n) if c._nodes[i].has_nondangling_edge()]
		m = len(a)
		if m > 6:
			return c
		if m < n:
			ia = [-1]*n
			for i in range(m):
				ia[a[i]] = i
			d = tc.Circuit(m)
			d.append(c, ia)
			c = tc.Circuit(n)
			c.any(*a, unitary=d.matrix())
			return c
		else:
			d = tc.Circuit(n)
			d.any(*range(n), unitary=c.matrix())
			return d
	return wrapper


def set_state(n: int, qubits: list) -> tc.Circuit:
    """
    Initialize a quantum circuit to represent an integer state.

    Parameters
    ----------
    n : int
        The integer value to encode in binary.
    qubits : list
        List of qubit indices ordered from most significant bit (MSB) to least significant bit (LSB).
        The binary representation of n is mapped onto these qubits.

    Returns
    -------
    tc.Circuit
        A circuit with the initial state |n> on the specified qubits.
    """
    c = tc.Circuit(max(qubits) + 1)
    l = 0
    while n > 0:
        if n % 2 == 1:
            c.X(qubits[-l - 1])
        else:
            c.I(qubits[-l - 1])
        n //= 2
        l += 1
    return c

@ block
def cset_state(p: int, n: int, qubits: list) -> tc.Circuit:
    """
    Controlled version of set_state, used in corner case of cond_ECC_base()

    Parameters
    ----------
    p : int
        The control qubit
    n : int
        The integer value to encode in binary.
    qubits : list
        List of qubit indices ordered from most significant bit (MSB) to least significant bit (LSB).
        The binary representation of n is mapped onto these qubits.

    Returns
    -------
    tc.Circuit
        A circuit with the initial state |n> on the specified qubits.
    """
    c = tc.Circuit(max(qubits + [p]) + 1)
    l = 0
    while n > 0:
        if n % 2 == 1:
            c.cnot(p, qubits[-l - 1])
        n //= 2
        l += 1
    return c


def out(state: np.ndarray, a=None, b=None) -> Union[int, list[int]]:
    """
    Extract sub-bitstring or non-zero amplitude indices from a quantum state.

    With one argument (state): returns a list (or single integer) of computational
    basis state indices whose amplitude magnitude > 0.001.

    With three arguments (state, a, b): extracts bits from position a to b
    (0-indexed, a <= b). For each basis state with non-zero amplitude, the substring
    is interpreted as a binary number, and the resulting integer(s) is/are returned.

    Parameters
    ----------
    state : np.ndarray
        State vector of the quantum circuit (length 2^n).
    a : int, optional
        Start bit position (0-indexed, more significant, must be <= b).
    b : int, optional
        End bit position (0-indexed, less significant, must be >= a).

    Returns
    -------
    int or list of int
        If a and b are None: returns the index (or list of indices) of basis states
        with |amplitude| > 0.001.
        If a and b are given: returns the integer formed by bits a..b,
        or a list of such integers if multiple basis states are present.
    """
    n_qubits = len(state).bit_length() - 1  # because state length = 2^n_qubits

    if a is None and b is None:
        # Return indices with significant amplitude (instead of printing)
        indices = [i for i, amp in enumerate(state) if abs(amp) > 0.001]
        if len(indices) == 1:
            return indices[0]
        else:
            return indices

    # Extract substring from bits a to b (0‑indexed, inclusive, a <= b)
    results = []
    for i, amp in enumerate(state):
        if abs(amp) > 0.001:
            binary = f"{i:0{n_qubits}b}"
            sub_bits = binary[a:b+1]   # b+1 because Python end is exclusive
            value = int(sub_bits, 2)
            results.append(value)

    if len(results) == 1:
        return results[0]
    else:
        return results

@ block
def ccphase(p: int, x: int, y: int, theta: float) -> tc.Circuit:
    """
    Controlled-controlled phase gate.

    Parameters
    ----------
    p : int
        First control qubit index.
    x : int
        Second control qubit index.
    y : int
        Target qubit index.
    theta : float
        Phase angle (in radians).

    Returns
    -------
    tc.Circuit
        A circuit implementing the CCPhase gate.
    """
    c = tc.Circuit(max([p] + [x] + [y]) + 1)
    c.cphase(x, y, theta=theta / 2.0)
    c.cnot(p, x)
    c.cphase(x, y, theta=-theta / 2.0)
    c.cnot(p, x)
    c.cphase(p, y, theta=theta / 2.0)
    return c

# Example usage:
# c = set_state(5, [0, 1, 2, 3])   # |0101> on qubits 0 (MSB) .. 3 (LSB)
# sv = c.state()
# out(sv)                           # prints indices with amplitude (e.g., 5)
# num = out(sv, 0, 3)              # extracts bits 0..3 = "0101" -> integer 5
# print(num)                        # outputs 5

