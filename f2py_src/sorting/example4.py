from f_sorting import sort
import numpy as np

sort.mb_breakpoint = 2
arr1 = np.random.rand(10)
arr2 = np.random.rand(10, 3)

print(arr2)
sort.sort_1d_dp_asc(arr2[:,1])
print(arr2)
arr2 = sort.apply_keytable_2d_dp(arr2)
print(arr2)





"""
arr2 = np.copy(arr)
print(arr)
arr_s = sort.sort_1d_dp_asc(arr)
print(arr_s)

print(sort.keytable)

print()
print(arr2)
arr2 = sort.apply_keytable_1d_dp(arr2)
print(arr2)
"""




