import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import tensorcircuit as tc
from shor_ecc.add_const import *

@ block
def mod_add(x: list, y: list, z: list) -> tc.Circuit:
    '''
    modular addition of two quantum registers.
    |x>|y> -> |x>|x+y(mod 7)>
    Here we assume x, y < 7.

    Parameters
    ----------
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    y : list
        Qubit indices for the 4-qubit register holding |y>.
    z : list
        Qubit indices for the ancilla register (at least 1 qubit).
    
    Returns
    -------
    tc.Circuit
        Circuit implementing the modular addition.
    '''
    c = tc.Circuit(max(x + y + z) + 1)
    c.append(add(x, y)) # x + y
    c.append(add_const(y, 9)) # x + y - 7
    c.cnot(y[0], z[0])
    c.append(cadd_const(z[0], y, 7)) # x + y if x + y < 7, x + y - 7 if x + y >= 7
    c.append(add(x, y).inverse()) # x if x + y < 7, x - 7 if x + y >= 7
    c.x(y[0])
    c.cnot(y[0], z[0])
    c.x(y[0]) # uncomputation of the ancilla
    c.append(add(x, y)) # x + y if x + y < 7, x + y - 7 if x + y >= 7
    return c



@ block
def cmod_add(p: int, x: list, y: list, z: list) -> tc.Circuit:
    '''
    Controlled modular addition.
    |p>|x>|y> -> |p>|x>|y + p*x (mod 7)>
    If control qubit p is |0>, nothing changes; if p is |1>, the modular addition is performed.
    Here we assume x, y < 7.

    Parameters
    ----------
    p : int
        Index of the control qubit.
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    y : list
        Qubit indices for the 4-qubit register holding |y>.
    z : list
        Qubit indices for the ancilla register (at least 1 qubit).
    
    Returns
    -------
    tc.Circuit
        Circuit implementing the controlled modular addition.
    '''
    c = tc.Circuit(max(x + y + z + [p]) + 1)
    c.append(add(x, y))
    c.append(cadd_const(p, y, 9))
    c.toffoli(p, y[0], z[0])
    c.append(ccadd_const(p, z[0], y, 7))
    c.append(add(x, y).inverse())
    c.x(y[0])
    c.toffoli(p, y[0], z[0])
    c.x(y[0])
    c.append(cadd(p, x, y))
    return c

