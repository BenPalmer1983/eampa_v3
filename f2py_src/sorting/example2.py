from f_sorting import sort
import numpy as np
import time

s = [10,50,100, 250, 1000]

for si in s:
  start = time.time()
  arr = np.random.rand(si)
  for i in range(10000):
    arr_n = np.copy(arr)
    arr_s = sort.sort_1d_dp_asc(arr_n)
  print(si, time.time() - start)
print(arr_s)

s = [10,50,100, 250, 1000]
for si in s:
  start = time.time()
  arr = np.random.randint(0,100,si)
  for i in range(10000):
    arr_n = np.copy(arr)
    arr_s = sort.sort_1d_int_asc(arr_n)
  print(si, time.time() - start)
print(arr_s)