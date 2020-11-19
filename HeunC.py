import cmath
from ctypes import *

HeunClib = cdll.LoadLibrary('HeunC_lib.so')
HeunClib.HeunC_func.argtypes = [c_double, c_double, c_double, c_double, c_double,\
                          c_double, c_double, c_double, c_double, c_double,\
                          c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double),\
                          POINTER(c_double), POINTER(c_int), POINTER(c_double)]

# Python class to compute the Heun function
class HeunC:
  def __init__(self,q,alpha,gamma,delta,epsilon):
    self.q = q
    self.alpha = alpha
    self.gamma = gamma
    self.delta = delta
    self.epsilon = epsilon
  def compute(self, z):
      # initialise values
      val_real = c_double(0.0)
      val_imag = c_double(0.0) 
      dval_real = c_double(0.0) 
      dval_imag = c_double(0.0) 
      numb_c = c_int(0)
      err_c = c_double(0.0)
      # Populate values
      HeunClib.HeunC_func(self.q.real,self.q.imag,self.alpha.real,self.alpha.imag,self.gamma.real,self.gamma.imag,self.delta.real,\
			self.delta.imag,self.epsilon.real,self.epsilon.imag,z,\
                          byref(val_real),byref(val_imag),byref(dval_real),byref(dval_imag),byref(numb_c),byref(err_c))
      val = complex(val_real.value,val_imag.value)
      dval = complex(dval_real.value,dval_imag.value)
      numb = numb_c.value
      err = err_c.value
      return (val, dval, numb, err)   
#   
