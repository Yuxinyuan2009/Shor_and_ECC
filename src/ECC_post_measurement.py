from ECC_base import *
import json
import sympy

def order(x: int, y: int, a: int) -> int:
    '''
    calculate the order of P(x, y). In the real-world sceneries, this order is given. 
    This order can be quickly calculated by Schoof-Elkies-Atkin (SEA) algorithm.

    Parameters
    ----------
    x: int
        x-coordinate of P
    y: int
        y-coordinate of Q
    a: int
        coefficient a in the elliptic curve equation
    
    Returns
    -------
    int
        the order of P(x, y)
    '''
    s, t = 0, 0
    for i in range(100):
        s, t = classical_ECC_add(x, y, s, t, a)
        if (s, t) == (0, 0):
            return i + 1

# print(order(2, 2, 0))

def get_ket_distribution(x1: int, y1: int, a: int) -> np.ndarray:
    with open("assets/output.json", 'r') as f:
        data = json.load(f)

    vec = np.array(data)

    r = order(x1, y1, a)

    s = np.zeros(r)
    for i in range(1024):
        m = i / 32
        n = i % 32
        x = round(m * r / 32)
        y = round(n * r / 32)
        s[y * (x ** (sympy.totient(r) - 1)) % r] += vec[i]
    
    return s

print(get_ket_distribution(2, 2, 0))