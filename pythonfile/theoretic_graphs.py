#!/usr/bin/env python3

from ctypes import cdll, c_double, CFUNCTYPE
import numpy as np
import seaborn
from pathlib import Path
import matplotlib
matplotlib.use("pgf")
pgf_with_pdflatex = {
    "pgf.texsystem": "pdflatex",
    "font.family": "serif", # use serif/main font for text elements
    "text.usetex": True,
    "errorbar.capsize": 0.5,
    "pgf.preamble": [
         r"\usepackage[utf8]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         r"\usepackage{mathpazo}",
             r"\usepackage[version-1-compatibility]{siunitx}"
         ]
}
matplotlib.rcParams.update(pgf_with_pdflatex)
from matplotlib import pyplot as plt
# import netCDF4

try:
  lib = cdll.LoadLibrary("pythonfile/libProjectPython.dylib")
except OSError as e:
  lib = cdll.LoadLibrary("pythonfile/libProjectPython.so")

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

def plot_contours(ufunc, title, filename):
  fig, plots = plt.subplots(2, 2, 'col', 'row',  True, gridspec_kw={
                              'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(5,5))
  t, d = np.meshgrid(temps, depths)
  cp = plots[0][0].contour(t, d, ufunc(d,t, 35))
  plt.clabel(cp, fmt="%d", rightside_up=False)
  plots[0][0].set_ylim(2000, 0)
  plots[0][0].set_ylabel("Depth (m)")
  plots[0][0].set_xlabel(r"Temperature (\si{\degreeCelsius})")
  
  s, d = np.meshgrid(salinities, depths)
  cp = plots[0][1].contour(s, d, ufunc(d, 10, s))
  plt.clabel(cp, fmt="%d", rightside_up=False)
  plots[0][1].set_ylim(2000, 0)
  plots[0][1].set_ylabel("Depth (m)")
  plots[0][1].set_xlabel("Salinity (ppt)")

  t,s = np.meshgrid(temps, salinities)
  cp = plots[1][0].contour(t,s, ufunc(1000, t, s))
  plt.clabel(cp, fmt="%d")
  plots[1][0].set_xlabel(r"Temperature (\si{\degreeCelsius})")
  plots[1][0].set_ylabel("Salinity (ppt)")
  fig.suptitle(title)

  plots[0][0].grid()
  plots[0][1].grid()
  plots[1][0].grid()
  fig.delaxes(plots[1][1])
  fig.savefig(filename)

plot_contours(ufunc_unesco, "UNESCO Equation (Chen and Millero 1995)", Path("final_output/figures/unesco.pdf"))
plot_contours(ufunc_leroy, "Leroy et al. 2008", Path("final_output/figures/leroy.pdf"))

plt.figure()
plt.plot(1000*ufunc_ideal(depths, 1160, 1.3, 1.45, 1.14e-3), depths)
plt.ylim(2000, 0)
plt.xlabel(r"Speed of Sound (\si{\meter\per\second})")
plt.ylabel(r"Depth (\si{\meter})")
plt.savefig(Path("final_output/figures/ideal.pdf"))