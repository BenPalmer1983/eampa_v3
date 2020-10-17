from f_sorting import bubblesort
from f_sorting import mergesort
from f_sorting import rng
from f_sorting import shuffle
from f_sorting import sort
import numpy as np
import time


arr = rng.randint_1d(1,50,1000000)
mergesort.set_1d_int(arr)
mergesort.ms_sort()
kt = mergesort.key_table
arr2 = mergesort.match_key_table_1d_int(arr)

print(mergesort.time)
print(arr2[0:20])
print(arr2[-10:-1])

arr = rng.randdp_1d(1.0,50.0,1000000)
mergesort.set_1d_dp(arr)
mergesort.ms_sort()
arr3 = mergesort.match_key_table_1d_dp(arr)


print(mergesort.time)
print(arr3[0:20])
print(arr3[-10:-1])
print(arr3[-1])


