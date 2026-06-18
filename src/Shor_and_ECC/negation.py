import tensorcircuit as tc
from Shor_and_ECC.mod_add import *

@ block
def negation(x: list, z: list) -> tc.Circuit:
    '''
    modular negation operation on a quantum register |x>.
    |x> -> |-x(mod 7)>
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
        Circuit implementing the modular negation operation.
    '''
    # Use the same method of mod_add
    c = tc.Circuit(max(x + z) + 1)
    c.x(x[0])
    c.x(x[1])
    c.x(x[2])
    c.x(x[3])
    c.append(add_const(x, 9))
    c.cnot(x[0], z[0])
    c.append(cadd_const(z[0], x, 9))
    c.append(add_const(x, 14))
    c.cnot(x[0], z[0])
    c.append(add_const(x, 1))
    return c



@ block
def cnegation(p: int, x: list, z: list) -> tc.Circuit:
    '''
    Controlled modular negation operation on a quantum register |x>.
    If control qubit p is 1, |x> -> |-x(mod 7)>, otherwise |x> remains unchanged.
    Here we assume x < 7.

    Parameters
    ----------
    p : int
        Qubit index for the control qubit.
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    z : list
        Qubit indices for the ancilla register (at least 1 qubit).

    Returns
    -------
    tc.Circuit
        Circuit implementing the controlled modular negation operation.
    '''
    c = tc.Circuit(max(x + z + [p]) + 1)
    c.cnot(p, x[0])
    c.cnot(p, x[1])
    c.cnot(p, x[2])
    c.cnot(p, x[3])
    c.append(cadd_const(p, x, 9))
    c.cnot(x[0], z[0])
    c.append(ccadd_const(p, z[0], x, 9))
    c.append(cadd_const(p, x, 14))
    c.toffoli(p, x[0], z[0])
    c.append(cadd_const(p, x, 1))
    return c

