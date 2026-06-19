import time
import tensorcircuit as tc
import numpy as np
from shor_ecc.ECC_base import *
import json

# Start timing
start_time = time.time()

x1, y1 = 2, 2
x2, y2 = 5, 4
a = 0

c = tc.Circuit(22)
for i in range(10):
    c.H(i)
for i in [4, 3, 2, 1, 0]:
    c.append(cond_ECC_add(i, x1, y1, [10, 11, 12, 13], [14, 15, 16, 17], a, [18, 19, 20, 21]))
    x1, y1 = classical_ECC_add(x1, y1, x1, y1, a)
for i in [9, 8, 7, 6, 5]:
    c.append(cond_ECC_add(i, x2, y2, [10, 11, 12, 13], [14, 15, 16, 17], a, [18, 19, 20, 21]))
    x2, y2 = classical_ECC_add(x2, y2, x2, y2, a)
c.append(QFT([0, 1, 2, 3, 4]).inverse())
c.append(QFT([5, 6, 7, 8, 9]).inverse())
print("yeah")
# about 110000 gates

vec = c.state()
def reduce_array(arr):
    arr = np.abs(arr)**2
    n = 22  # total number of qubits
    m = n - 10  # number of qubits to sum over
    chunk_size = 2 ** m  # size of each block
    # Reshape to (2^10, 2^m) and sum over the second axis
    reshaped = arr.reshape(-1, chunk_size)
    reduced = np.sum(reshaped, axis=1)
    return reduced

vec = reduce_array(vec)

with open("assets/output.json", "w") as f:
    json.dump(vec.tolist(), f)

# End timing and display runtime
end_time = time.time()
print(f"Runtime: {end_time - start_time:.4f} seconds")
# Runtime: 197 seconds(for 5 * 2 = 10 qubits in first register)
# Runtime: 90 seconds(for 4 * 2 = 8 qubits in first register)