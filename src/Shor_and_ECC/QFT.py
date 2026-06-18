import tensorcircuit as tc
import numpy as np
from Shor_and_ECC.base import *

@ block
def QFT(qubits) -> tc.Circuit:
    """
    Construct a Quantum Fourier Transform circuit on a given set of qubits.

    Parameters
    ----------
    qubits : list of int
        List of qubit indices on which to apply the QFT.
        The order matters: qubits[0] is the most significant qubit (|x1>)
        and qubits[-1] is the least significant (|xm>).

    Returns
    -------
    tc.Circuit
        A TensorCircuit object representing the QFT gate.
    """
    # Number of qubits needed is the maximum index + 1 (in case indices are not contiguous)
    n = max(qubits) + 1
    m = len(qubits)          # number of qubits in the QFT
    c = tc.Circuit(n)        # create empty circuit with enough qubits

    # Main QFT loop: iterate over each qubit from most significant to least
    for i in range(m):
        c.H(qubits[i])
        for j in range(i + 1, m):
            c.cphase(qubits[j], qubits[i], theta=np.pi / (2 ** (j - i)))

    # Bit-reversal permutation: swap the i-th and (m-1-i)-th qubits
    # This corrects the order of the output amplitudes.
    for i in range(m // 2):
        c.swap(qubits[i], qubits[m - 1 - i])

    return c

# Example usage:
# qft_circuit = QFT([0, 1, 2])   # 3-qubit QFT on qubits 0,1,2