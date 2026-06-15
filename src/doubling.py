import tensorcircuit as tc
from mod_add import *

@ block
def doubling(x: list, z: list) -> tc.Circuit:
    '''
    Doubling operation on a quantum register |x>.
    |x> -> |2x(mod p)>
    Here we assume x < 7.

    Parameters
    ----------
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    z : list
        Qubit indices for the ancilla register (at least 1 qubit).
    
    Returns
    -------
    tc.Circuit
        Circuit implementing the doubling operation.
    '''
    c = tc.Circuit(max(x + z) + 1)
    for i in range(3):
        c.swap(x[i], x[i + 1])
    # |x> -> |2x>
    c.cnot(x[0], z[0])
    c.append(cadd_const(z[0], x, 9))
    c.cnot(x[3], z[0])
    return c
