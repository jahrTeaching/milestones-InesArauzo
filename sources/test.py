from numpy import array
from resources.Math_Operators import LU_Inverse

A = array ([[1,2], [3,4]])
invA = LU_Inverse (A)
print(invA)