
from sorting import sorting
from f_sorting import bubblesort
from f_sorting import mergesort
from f_sorting import rng



arr1 = rng.randint_1d(1,50,100)
sorting.sort(arr1)


arr1 = rng.randdp_1d(1,50,100)
sorting.sort(arr1)


"""
arr1 = rng.randint_1d(1,50,100)
arr2 = rng.randint_1d(1,50,100)

print(arr1)
print(arr2)

mergesort.set_1d_int(arr1)
mergesort.ms_sort()
kt = mergesort.key_table
arr1 = mergesort.match_key_table_1d_int(arr1)
arr2 = mergesort.match_key_table_1d_int(arr2)

print(arr1)
print(arr2)
"""