import tensorcircuit as tc
from shor_ecc.mod_add import *
from shor_ecc.doubling import *
from shor_ecc.mod_add_const import *
from tqdm import tqdm

@ block
def mod_mul(x: list, y: list, z: list) -> tc.Circuit:
    '''
    modular multiplication operation on quantum registers |x> and |y>.
    |x>|y>|z>  ->|x>|y>|z+x*y(mod 7)>
    Here we assume x, y, z < 7.

    Parameters
    ----------
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    y : list
        Qubit indices for the 4-qubit register holding |y>.
    z : list
        Qubit indices for the 4-qubit register holding |z>.
    
    Returns
    -------
    tc.Circuit
        Circuit implementing the modular multiplication operation.
    '''
    c = tc.Circuit(max(x + y + z) + 1)
    for i in range(3):
        c.append(cmod_add(x[3 - i], y, z, [x[0]])) # add x[i] * y to z
        c.append(doubling(y, [x[0]]))
    return c

# Test
'''
for i in tqdm(np.arange(7)):
    for j in range(7):
        c = tc.Circuit(12)
        c.append(set_state(i, [0, 1, 2, 3]))
        c.append(set_state(j, [4, 5, 6, 7]))
        c.append(mod_mul([0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]))
        assert out(c.state(), 8, 11) == (i * j) % 7, f"test mod_mul with parameters {i}, {j} has problem, expected {(i * j) % 7}, got {out(c.state(), 8, 11)}"
'''

@ block
def mod_square(x: list, y: list) -> tc.Circuit:
    '''
    modular squaring operation on a quantum register |x>.
    |x>|y> -> |x>|y+x^2(mod 7)>
    Here we assume x, y < 7.

    Parameters
    ----------
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    y : list
        Qubit indices for the 4-qubit register holding |y>.
    
    Returns
    -------
    tc.Circuit
        Circuit implementing the modular squaring operation.
    '''
    # Note that this realization is O(n^2). We apply this method to minimize ancilla qubtis usage.
    # There is a more efficient realization using the Karatsuba algorithm, which is left for future work.
    c = tc.Circuit(max(x + y) + 1)
    for i in range(1, 4):
        for j in range(i, 4):
            if i == j:
                c.append(cmod_add_const(x[i], y, [x[0]], 2**(6 - i - j) % 7)) # add x[i] * x[i]
            else:
                c.append(ccmod_add_const(x[i], x[j], y, [x[0]], 2**(7 - i - j) % 7)) # add 2 * x[i] * x[j]
    return c



@ block
def cmod_square(p: int, x: list, y: list) -> tc.Circuit:
    '''
    controlled modular squaring operation on a quantum register |x>.
    |x>|y> -> |x>|y+x^2(mod p)>
    If control qubit p is |0>, nothing changes; if p is |1>, the modular squaring is performed.
    Here we assume x, y < 7.

    Parameters
    ----------
    p : int
        Index of the control qubit.
    x : list
        Qubit indices for the 4-qubit register holding |x>.
    y : list
        Qubit indices for the 4-qubit register holding |y>.
    
    Returns
    -------
    tc.Circuit
        Circuit implementing the controlled modular squaring operation.
    '''
    c = tc.Circuit(max([p] + x + y) + 1)
    for i in range(1, 4):
        for j in range(i, 4):
            if i == j:
                c.append(ccmod_add_const(p, x[i], y, [x[0]], 2**(6 - i - j) % 7))
            else:
                c.append(cccmod_add_const(p, x[i], x[j], y, [x[0]], 2**(7 - i - j) % 7))
    return c



@ block
def inv(x):
    c = tc.Circuit(max(x) + 1)
    c.swap(x[1], x[2])
    return c

