import numpy as np
import cmath
import ctypes

HeunClib = ctypes.cdll.LoadLibrary('HeunC_lib.so')
Kerrlib.Rfunc.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,\
                          ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,\
                          ctypes.c_double, POINTER(ctypes.c_double), POINTER(ctypes.c_double), POINTER(ctypes.c_double),\
                          POINTER(ctypes.c_double), POINTER(ctypes.c_int), POINTER(ctypes.c_double)]

class HeunC(q,alpha,gamma,delta,epsilon):
  def __init__(self,q,alpha,gamma,delta,epsilon):
    self.q = q
    self.alpha = alpha
    self.gamma = gamma
    self.delta = delta
    self.epsilon = epsilon
  def compute_single(self, z):
      # initialise values
      val_real = ctypes.c_double(0.0)
      val_imag = ctypes.c_double(0.0) 
      dval_real = ctypes.c_double(0.0) 
      dval_imag = ctypes.c_double(0.0) 
      numb = ctypes.c_int(0)
      err = ctypes.c_double(0.0)
      # Populate values
      HeunClib.
