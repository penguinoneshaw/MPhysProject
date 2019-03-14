#!/usr/bin/env python3

from ctypes import cdll, c_double, CFUNCTYPE
import numpy as np
import seaborn
from matplotlib import pyplot as plt
# import netCDF4

try:
  lib = cdll.LoadLibrary("./libProjectPython.dylib")
except OSError as e:
  lib = cdll.LoadLibrary("./libProjectPython.so")

lib.unesco_depth.argtypes = [c_double,c_double,c_double]
lib.unesco_depth.restype = c_double

lib.unesco_pressure.argtypes = [c_double,c_double,c_double]
lib.unesco_pressure.restype = c_double

lib.leroy_et_al.argtypes = [c_double,c_double,c_double]
lib.leroy_et_al.restype = c_double

lib.ideal_sound_channel.argtypes = [c_double,c_double,c_double,c_double,c_double]
lib.ideal_sound_channel.restype = c_double

depths = np.linspace(0,2000,5000,dtype=np.double)
pressures = np.linspace(0,1000,5000,dtype=np.double)
temps = np.linspace(0, 40, 100, dtype=np.double)
salinities = np.linspace(0, 40, 100, dtype=np.double)

ufunc_unesco = np.frompyfunc(lib.unesco_depth, 3, 1)
ufunc_leroy = np.frompyfunc(lib.leroy_et_al, 3, 1)
ufunc_ideal = np.frompyfunc(lib.ideal_sound_channel, 5, 1)

def plot_contours(ufunc, title):
  plt.figure()
  plt.subplot(221)
  t, d = np.meshgrid(temps, depths)
  cp = plt.contour(t, d, ufunc(d,t, 35))
  plt.clabel(cp, fmt="%d", rightside_up=False)
  plt.ylim(2000, 0)
  plt.ylabel("Depth (m)")
  plt.xlabel("Temperature (degrees C)")

  plt.subplot(222)
  s, d = np.meshgrid(salinities, depths)
  cp = plt.contour(s, d, ufunc(d, 10, s))
  plt.clabel(cp, fmt="%d", rightside_up=False)
  plt.ylim(2000, 0)
  plt.ylabel("Depth (m)")
  plt.xlabel("Salinity (ppt)")

  plt.subplot(223)
  t,s = np.meshgrid(temps, salinities)
  cp = plt.contour(t,s, ufunc(1000, t, s))
  plt.clabel(cp, fmt="%d")
  plt.xlabel("Temp (degrees Celsius)")
  plt.ylabel("Salinity (ppt)")
  plt.suptitle(title)

plot_contours(ufunc_unesco, "UNESCO Equation (Chen and Millero 1998)")
plot_contours(ufunc_leroy, "Leroy et. al. 2008")

plt.figure()
plt.plot(1000*ufunc_ideal(depths, 1160, 1.3, 1.45, 1.14e-3), depths)
plt.ylim(2000, 0)
plt.xlabel("Speed of Sound (m/s)")
plt.ylabel("Depth (m)")
plt.show()