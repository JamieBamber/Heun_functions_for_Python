import numpy as np
import cmath
import HeunC

# Initialise HeunC class
HC = HeunC.HeunC(0, -0.5, -0.5, 0, 1.0/8)
"""
With these parameter choices HeunC(0, -0.5, -0.5, 0, 1.0/8, 1-x) = sqrt(x)
Note that the HeunC branch cut means we can't have real x values less than 0
"""

def test_HeunC_function(r):
  output = np.zeros(len(r))
  for i in range(0, len(r)):
      val, dval, numb, err = HC.compute(1-r)
      output[i] = val.real
  return output

r = [1, 4, 9, 16, 25]
print(test_HeunC_function(r))
