import numpy as np
import cmath
import ctypes

HeunClib = ctypes.cdll.LoadLibrary('HeunC_lib.so')
Kerrlib.Rfunc.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,\
                          ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,\
                          ctypes.c_double, POINTER(ctypes.c_double), POINTER(ctypes.c_double), POINTER(ctypes.c_double),\
                          POINTER(ctypes.c_double), POINTER(ctypes.c_int), POINTER(ctypes.c_double)]

# Python class to compute the Heun function
class HeunC(q,alpha,gamma,delta,epsilon):
  def __init__(self,q,alpha,gamma,delta,epsilon):
    self.q = q
    self.alpha = alpha
    self.gamma = gamma
    self.delta = delta
    self.epsilon = epsilon
  def compute(self, z):
      # initialise values
      val_real = ctypes.c_double(0.0)
      val_imag = ctypes.c_double(0.0) 
      dval_real = ctypes.c_double(0.0) 
      dval_imag = ctypes.c_double(0.0) 
      numb_c = ctypes.c_int(0)
      err_c = ctypes.c_double(0.0)
      # Populate values
      HeunClib.HeunC_func(q.real,q.imag,alpha.real,alpha.imag,gamma.real,gamma.imag,delta.real,delta.imag,epsilon.real,epsilon.imag,z,\
                          ctypes.byref(val_real),ctypes.byref(val_imag),ctypes.byref(dval_real),ctypes.byref(dval_imag),ctypes.byref(numb_c),ctypes.byref(err_c))
      val = complex(val_real.value,val_imag.value)
      dval = complex(dval_real.value,dval_imag.value)
      numb = numb_c.value
      err = err_c.value
      return (val, dval, numb, err)   
#   
