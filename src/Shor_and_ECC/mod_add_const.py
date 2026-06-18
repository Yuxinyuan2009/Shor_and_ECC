import tensorcircuit as tc
from Shor_and_ECC.add_const import *

@ block
def mod_add_const(x: list, z: list, t: int) -> tc.Circuit:
    '''
    modular addition of a constant t to a quantum register |x>.
    |x> -> |x+t(mod 7)>
    Here we assume x < 7.

    Parameters
    ----------
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    z : list
        Qubit indices for the ancilla register (at least 1 qubit).
    t : int
        Constant to add modulo 7.
    
    Returns
    -------
    tc.Circuit
        Circuit implementing the modular addition.
    '''
    t = t % 7
    c = tc.Circuit(max(x + z) + 1)
    c.append(add_const(x, t + 1))
    c.cnot(x[0], z[0])
    c.append(cadd_const(z[0], x, 9))
    c.append(add_const(x, -(t + 1)))
    c.cnot(x[0], z[0])
    c.append(add_const(x, t))
    return c


@ block
def cmod_add_const(p: int, x: list, z: list, t: int) -> tc.Circuit:
    '''
    Controlled modular addition of a constant t.
    |p>|x> -> |p>|x + p*t (mod 7)>
    If control qubit p is |0>, nothing changes; if p is |1>, the modular addition is performed.
    Here we assume x < 7. 

    Parameters
    ----------
    p : int
        Index of the control qubit.
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    z : list
        Qubit indices for the ancilla register (at least 1 qubit).
    t : int
        Constant to add modulo 7.

    Returns
    -------
    tc.Circuit
        Circuit implementing the controlled modular addition.
    '''
    t = t % 7
    c = tc.Circuit(max(x + z + [p]) + 1)
    c.append(add_const(x, t + 1))
    c.toffoli(p, x[0], z[0])
    c.append(ccadd_const(p, z[0], x, 9))
    c.append(add_const(x, -(t + 1)))
    c.toffoli(p, x[0], z[0])
    c.append(cadd_const(p, x, t))
    return c



@ block
def ccmod_add_const(p: int, q: int, x: list, z: list, t: int) -> tc.Circuit:
    '''
    Controlled-controlled modular addition of a constant t.
    |p>|q>|x> -> |p>|q>|x + p*q*t (mod 7)>
    If both control qubits p and q are |1>, the modular addition is performed.
    Here we assume x < 7.

    Parameters
    ----------
    p : int
        Index of the first control qubit.
    q : int
        Index of the second control qubit.
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    z : list
        Qubit indices for the ancilla register (at least 1 qubit).
    t : int
        Constant to add modulo 7.

    Returns
    -------
    tc.Circuit
        Circuit implementing the controlled-controlled modular addition.
    '''
    # Here we use the factorization of multicontrol gates.
    t = t % 7
    c = tc.Circuit(max(x + z + [p] + [q]) + 1)
    c.append(cmod_add_const(p, x, z, t * 4 % 7))
    c.cnot(q, p)
    c.append(cmod_add_const(p, x, z, t * 3 % 7))
    c.cnot(q, p)
    c.append(cmod_add_const(q, x, z, t * 4 % 7))
    return c



@ block
def cccmod_add_const(p: int, q: int, r: int, x: list, z: list, t: int) -> tc.Circuit:
    '''
    Controlled-controlled-controlled modular addition.
    |p>|q>|r>|x> -> |p>|q>|r>|x + p*q*r*x (mod 7)>
    If control qubits p, q, r are all |1>, the modular addition is performed.
    Here we assume x < 7.

    Parameters
    ----------
    p : int
        Index of the first control qubit.
    q : int
        Index of the second control qubit.
    r : int
        Index of the third control qubit.
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    z : list
        Qubit indices for the ancilla register (at least 1 qubit).
    t : int
        Constant to add modulo 7.

    Returns
    -------
    tc.Circuit
        Circuit implementing the controlled-controlled-controlled modular addition.
    '''
    t = t % 7
    c = tc.Circuit(max([p, q, r] + x + z) + 1)
    c.append(ccmod_add_const(p, r, x, z, t * 4 % 7))
    c.cnot(q, p)
    c.append(ccmod_add_const(p, r, x, z, t * 3 % 7))
    c.cnot(q, p)
    c.append(ccmod_add_const(q, r, x, z, t * 4 % 7))
    return c

