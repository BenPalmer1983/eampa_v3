from f_rng import rng
import numpy as np




bin = np.zeros(10)

for i in range(1000000):
  x = rng.randint(10,19)
  bin[x-10] = bin[x-10] + 1
print(bin)



#r2d = np.zeros((1000,10), dtype=np.int32)
r2d = rng.randint_2d(10, 20, 100, 1000)

print(r2d[1:4,1:4])
