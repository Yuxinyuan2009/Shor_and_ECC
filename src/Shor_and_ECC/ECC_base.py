import tensorcircuit as tc
from Shor_and_ECC.mod_mul import *
from Shor_and_ECC.negation import *

def classical_ECC_add(x1: int, y1: int, x2: int, y2: int, a: int) -> tuple[int, int]:
    '''
    classical implementation of point addition on elliptic curve defined by y^2 = x^3 + ax + b (mod 7).
    Two points P = (x1, y1), Q = (x2, y2) are given, return P + Q.
    Here P != O.
    Point addition followed Rules:
    - If Q != O:
        lambda = (y2 - y1) / (x2 - x1) (mod 7) if P != Q,

        = (3 * x1^2 + a) / (2 * y1) (mod 7) if P = Q

        x3 = lambda^2 - x1 - x2 (mod 7)
        y3 = lambda * (x1 - x3) - y1 (mod 7)
    - If Q = O:
        Q + P = P
    - If Q = -P:
        Q + P = O
    - If Q = -2P:
        Q + P = -P
    
    Parameters
    ----------
    x1: int
        x coordinate of P
    y1: int
        y coordinate of P
    x2: list
        x coordinate of Q, represented as a list of 4 qubits (little endian)
    y2: list
        y coordinate of Q, represented as a list of 4 qubits (little endian)
    a: int
        coefficient a in the elliptic curve equation
    
    Returns
    -------
    tuple[int, int]
        the coordinates of P + Q
    '''
    # Corner Case: Q = O, Q = -P, Q = P, Q = -2P
    if (x2, y2) == (0, 0):
        return x1, y1
    if (x2, y2) == (x1, (7 - y1) % 7):
        return 0, 0
    # Now we compute 2P
    l = (3 * x1**2 + a) * pow(2 * y1, 5, 7) % 7
    x3 = (l**2 - 2 * x1) % 7
    y3 = (l * (x1 - x3) - y1) % 7

    if (x2, y2) == (x1, y1):
        return x3, y3
    if (x2, y2) == (x3, (7 - y3) % 7):
        return x1, -y1 % 7
    # Now we compute P + Q
    l = (y2 - y1) * pow(x2 - x1, 5, 7) % 7
    x3 = (l**2 - x1 - x2) % 7
    y3 = (l * (x1 - x3) - y1) % 7
    return x3, y3

    
def ECC_add(x1: int, y1: int, x2: list, y2: list, a: int, z: list) -> tc.Circuit:
    '''
    Implementation of point addition on elliptic curve defined by y^2 = x^3 + ax + b (mod 7).
    P = (x1, y1), Q = (x2, y2)
    |Q> -> |Q + P>
    We use the following formulas for point addition:
    lambda = (y2 - y1) / (x2 - x1) (mod 7) if P != Q
    lambda = (3 * x1^2 + a) / (2 * y1) (mod 7) if P = Q
    x3 = lambda^2 - x1 - x2 (mod 7)
    y3 = lambda * (x1 - x3) - y1 (mod 7)
    Note: in this implememtation we don't consider corner case: Q = O, Q = -P, Q = -2P. We will handle them in cond_ECC_add.
    '''
    c = tc.Circuit(max(x2 + y2 + z) + 1)
    # We consider the state of |Q>, starting with: x2, y2, 0

    c.append(mod_add_const(x2, [z[0]], (7 - x1) % 7))
    c.append(mod_add_const(y2, [z[0]], (7 - y1) % 7))
    # now: x2 - x1, y2 - y1, 0
    c.multicontrol(*(x2[1:4] + y2[1:4] + [x2[0]]), unitary=tc.gates._x_matrix, ctrl=[0, 0, 0, 0, 0, 0])
    c.append(cmod_add_const(x2[0], z, [y2[0]], (3 * x1 ** 2 + a) * ((2 * y1)**5) % 7))
    c.multicontrol(*(x2[1:4] + y2[1:4] + [x2[0]]), unitary=tc.gates._x_matrix, ctrl=[0, 0, 0, 0, 0, 0])
    # corner case: P = Q
    c.append(inv(x2))
    c.append(mod_mul(x2, y2, z))
    c.append(inv(x2))
    # now: x2 - x1, y2 - y1, lambda
    c.append(mod_add_const(x2, [z[0]], x1 % 7))
    c.append(negation(x2, [z[0]]))
    c.append(mod_mul(x2, z, y2))
    c.append(negation(y2, [z[0]]))
    # now: -x2, lambda * x2 - y2 + y1, lambda
    c.append(mod_square(z, x2))
    c.append(mod_add_const(x2, [z[0]], (7 - x1) % 7))
    # now: x3, lambda * x2 - y2 + y1, lambda
    # = x3, lambda * x1, lambda
    c.append(negation(z, [x2[0]]))
    c.append(mod_mul(x2, z, y2))
    c.append(mod_add_const(x2, [z[0]], (7 - x1) % 7))
    c.append(negation(z, [x2[0]]))
    # now: x3 - x1, lambda * (x1 - x3), lambda
    c.append(inv(x2))
    c.append(mod_mul(x2, y2, z))
    c.append(inv(x2))
    # now: x3 - x1, lambda * (x1 - x3), 0
    c.append(mod_add_const(x2, [z[0]], x1 % 7))
    c.append(mod_add_const(y2, [z[0]], (7 - y1) % 7))
    return c

# Test
'''
c = tc.Circuit(12)
c.append(set_state(3, [0, 1, 2, 3]))
c.append(set_state(3, [4, 5, 6, 7]))
c.append(ECC_add(5, 6, [0, 1, 2, 3], [4, 5, 6, 7], 2, [8, 9, 10, 11]))
assert (out(c.state(), 0, 3), out(c.state(), 4, 7)) == classical_ECC_add(3, 3, 5, 6, 2)
'''

def cond_ECC_add(p: int, x1: int, y1: int, x2: list, y2: list, a: int, z: list) -> tc.Circuit:
    '''
    Implementation of controlled point addition on elliptic curve defined by y^2 = x^3 + ax + b (mod 7).
    P = (x1, y1), Q = (x2, y2), P != O
    |Q> -> |Q + P>
    We use the following formulas for point addition:
    - If Q != O:
        lambda = (y2 - y1) / (x2 - x1) (mod 7) if P != Q,

        = (3 * x1^2 + a) / (2 * y1) (mod 7) if P = Q

        x3 = lambda^2 - x1 - x2 (mod 7)
        y3 = lambda * (x1 - x3) - y1 (mod 7)
    - If Q = O:
        Q + P = P
    - If Q = -P:
        Q + P = O
    - If Q = -2P:
        Q + P = -P
    
    Parameters
    ----------
    p: int
        control qubit
    x1: int
        x coordinate of P
    y1: int
        y coordinate of P
    x2: list
        x coordinate of Q, represented as a list of 4 qubits (little endian)
    y2: list
        y coordinate of Q, represented as a list of 4 qubits (little endian)
    a: int
        coefficient a in the elliptic curve equation
    z: list
        auxiliary qubits, represented as a list of 4 qubits
    
    Returns
    -------
    tc.Circuit
        the circuit implementing the controlled point addition
    '''
    c = tc.Circuit(max([p] + x2 + y2 + z) + 1)
    # We consider the state of |Q>, starting with: x2, y2, 0

    c.append(mod_add_const(x2, [z[0]], (7 - x1) % 7))
    c.append(mod_add_const(y2, [z[0]], (7 - y1) % 7))
    # now: p = 1: x2 - x1, y2 - y1, 0
    # p = 0: x2 - x1, y2 - y1, 0
    c.multicontrol(*([p] + x2[1:4] + y2[1:4] + [x2[0]]), unitary=tc.gates._x_matrix, ctrl=[1, 0, 0, 0, 0, 0, 0])
    c.append(cmod_add_const(x2[0], z, [y2[0]], (3 * x1 ** 2 + a) * ((2 * y1)**5) % 7))
    c.multicontrol(*([p] + x2[1:4] + y2[1:4] + [x2[0]]), unitary=tc.gates._x_matrix, ctrl=[1, 0, 0, 0, 0, 0, 0])
    # corner case: P = Q
    c.append(inv(x2))
    c.append(mod_mul(x2, y2, z))
    c.append(inv(x2))
    # now: p = 1: x2 - x1, y2 - y1, lambda
    # p = 0: x2 - x1, y2 - y1, lambda
    c.append(mod_add_const(x2, [z[0]], x1 % 7))
    c.append(negation(x2, [z[0]]))
    c.append(mod_mul(x2, z, y2))
    c.append(negation(y2, [z[0]]))
    # now: p = 1: -x2, lambda * x2 - y2 + y1, lambda
    # p = 1 & Q = -P: -x1, 2y1, 0
    # p = 1 & Q = O: 0, y1, lambda, lambda = y1 / x1
    # p = 0: -x1, lambda * x2 - y2 + y1, lambda
    c.append(cmod_square(p, z, x2))
    c.append(cmod_add_const(p, x2, [z[0]], (7 - x1) % 7))
    # now: p = 1: x3, lambda * x2 - y2 + y1, lambda
    # = x3, lambda * x1, lambda
    # p = 1 & Q = -P: -2x1, 2y1, 0
    # p = 1 & Q = O: lambda^2 - x1, y1, lambda
    # p = 0: -x2, lambda * x2 - y2 + y1, lambda
    c.append(cnegation(p, z, [x2[0]]))
    c.append(mod_mul(x2, z, y2))
    c.x(p)
    c.append(cnegation(p, x2, [z[0]]))
    c.x(p)
    c.append(mod_add_const(x2, [z[0]], (7 - x1) % 7))
    c.append(cnegation(p, z, [x2[0]]))
    # now: p = 1: x3 - x1, lambda * (x1 - x3), lambda
    # p = 1 & Q = -P: -3x1, 2y1, 0
    # p = 1 & Q = O: lambda^2 - 2x1, -lambda * (lambda^2 - x1) + y1, lambda
    # p = 0: x2 - x1, -y2 + y1, lambda
    c.multicontrol(*([p] + x2[1:4] + [y2[0]]), unitary=tc.gates._x_matrix, ctrl=[1, 0, 0, 0])
    c.append(cmod_add_const(y2[0], z, [x2[0]], (7 - ((3 * x1**2 + a) * (2 * y1)**5 % 7)) % 7))
    c.multicontrol(*([p] + x2[1:4] + [y2[0]]), unitary=tc.gates._x_matrix, ctrl=[1, 0, 0, 0])
    # corner case: Q = -2P
    # 0, 0, lambda -> 0, 0, 0
    c.append(mod_add_const(y2, [z[0]], (7 - 2 * y1 % 7) % 7))
    c.multicontrol(*([p] + z[1:4] + y2[1:4] + [z[0]]), unitary=tc.gates._x_matrix, ctrl=[1, 0, 0, 0, 0, 0, 0])  # Q = -P: -3x1, 0, 0
    c.append(cset_state(z[0], -3 * x1 % 7, x2))                                                                 # Q = -P: 0, 0, 0
    c.append(cset_state(z[0], -x1 % 7, x2))                                                                     # Q = -P: -x1, 0, 0
    c.append(cset_state(z[0], -y1 % 7, y2))                                                                     # Q = -P: -x1, -y1, 0
    c.append(cset_state(z[0], y1 * (x1 ** 5) % 7, z[1:4]))                                                      # Q = -P: -x1, -y1, y1 / x1
    c.append(mod_add_const(y2, [x2[0]], y1 % 7))                                                                # change y2 to 0, so that z[0](x2[0]) can return to 0
    c.swap(z[0], y2[0])                                                                                         # swap the ancilla qubit to avoid affecting z(in the next step)
    c.append(mod_add_const(z, [x2[0]], 7 - y1 * (x1**5) % 7))                                                   # change y2 to 0, so that z[0](x2[0]) can return to 0
    c.multicontrol(*([p] + z[1:4] + y2[1:4] + [y2[0]]), unitary=tc.gates._x_matrix, ctrl=[1, 0, 0, 0, 0, 0, 0]) # reset x2[0](z[0])
    c.append(mod_add_const(y2, [x2[0]], y1 % 7))
    c.append(mod_add_const(z, [x2[0]], y1 * (x1**5) % 7))
    # corner case: Q = -P
    # -3x1, 2y1, 0 -> -x1, y1, y1 / x1
    c.append(inv(x2))
    c.append(mod_mul(x2, y2, z))
    c.append(inv(x2))
    # now: p = 1: x3 - x1, lambda * (x1 - x3), 0
    # p = 1 & Q = -P: -x1, y1, 0
    # p = 1 & Q = O: lambda^2 - 2x1, -lambda * (lambda^2 - 2x1), 0
    # p = 0: x2 - x1, -y2 + y1, 0
    c.append(mod_add_const(x2, [z[0]], x1 % 7))
    c.append(mod_add_const(y2, [z[0]], (7 - y1) % 7))
    c.x(p)
    c.append(cnegation(p, y2, [z[0]]))
    c.x(p)
    # now: p = 1: x3, y3, 0
    # p = 1 & Q = -P: 0, 0, 0
    # p = 1 & Q = O: lambda^2 - x1, lambda * (lambda^2 - x1), 0
    # p = 0: x2, y2, 0

    l = y1 * (x1**5) % 7    # lambda
    c.append(mod_add_const(x2, [z[0]], (x1 - l**2 % 7) % 7))
    c.append(mod_add_const(y2, [z[0]], (l**3 - l * x1) % 7))
    c.multicontrol(*([p] + x2[1:4] + y2[1:4] + [z[0]]), unitary=tc.gates._x_matrix, ctrl=[1, 0, 0, 0, 0, 0, 0]) # Q = O: 0, 0, 0
    c.append(cset_state(z[0], (2 * x1 - l**2) % 7, x2))
    c.append(cset_state(z[0], (l**3 - l * x1 + y1) % 7, y2))
    c.append(mod_add_const(x2, [y2[0]], (l**2 - 2 * x1) % 7))
    c.append(mod_add_const(y2, [x2[0]], (l * x1 - l**3 - y1) % 7))
    c.multicontrol(*([p] + x2[1:4] + y2[1:4] + [z[0]]), unitary=tc.gates._x_matrix, ctrl=[1, 0, 0, 0, 0, 0, 0])
    c.append(mod_add_const(x2, [z[0]], x1))
    c.append(mod_add_const(y2, [z[0]], y1))
    # corner case: Q = O
    # lambda^2 - x1, lambda * (lambda^2 - x1), 0 -> x1, y1, 0
    return c

# 13 qubits
# ~ 11000 gates


