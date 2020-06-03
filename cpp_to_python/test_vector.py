import test_vector
import numpy as np

a = np.zeros(10,np.float)
b = np.ones(10,np.float)
c = test_vector.vec_add(a,b)

a = [1,2,3]
b = [3,2,1]
c = test_vector.vec_add(a,b)

