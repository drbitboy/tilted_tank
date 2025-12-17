"""
Calculate volume of liquid in a flat-bottomed cylindrical tank resting
on a sloped base so the tank's cylindrical axis is tilted at angle phi
from vertical and the flat bottom is at angle phi from horizontal.

"""
import os
import sys
import numpy

class TILTED_TANK:
  def __init__(self, R=1.0, phi=30, n=64, **kwargs):
    """
  R - radius of tank
phi - angle of tilt, degrees
  n - number of voxels across the diameter of the tank
"""
    self.R, self.phideg, self.n = float(R), float(phi), int(n)
    self.phi = self.phideg * numpy.pi / 180.0
    self.cos,self.sin = numpy.cos(self.phi), numpy.sin(self.phi)
    self.tan = numpy.tan(self.phi)
    self.upperlim = 2.0 * self.R * self.sin
    self.halfedge = self.R / self.n
    self.edge = 2.0 * self.halfedge
    self.rowxs = (self.edge * numpy.arange(self.n)) + self.halfedge
    ylens = 2.0 * numpy.sqrt((self.R**2) - ((self.R-self.rowxs)**2))
    self.rowvols = (self.edge**2) * ylens
    self.rowbaseys = self.rowxs * self.sin
    self.rowdy = self.edge * self.cos

  def __repr__(self):
    return f"<{self.__class__.__name__} R:{self.R} phi:{self.phideg}deg n:{self.n}>"

  def volest(self, Level):

    Llim = max([0.0,min([Level,self.upperlim])])
    Vest,L = 0.0,Llim - (self.rowdy / 2.0)

    while L > 0.0:
      iw = numpy.where(L > self.rowbaseys)
      if not len(iw[0]): break
      Vest += numpy.sum(self.rowvols[iw])
      L -= self.rowdy

    return Vest

  def volcalc(self, Level):

    Llim = max([0.0,min([Level,self.upperlim])])

    H = Llim / self.sin
    X = (self.R - H) / self.R
    if X   >  1.0: X,H =  1.0,0.0
    elif X < -1.0: X,H = -1.0,2.0 * self.R
    Y = (2.0 * self.R * H) - (H**2)

    Xterm = (self.R**3) * (numpy.sqrt(1.0 - (X**2)) - (X*numpy.arccos(X)))
    Yterm = - (Y**1.5) / 3.0
    Vcalc = self.tan * (Xterm + Yterm)

    return Vcalc

module_kwargs = dict()
for arg in sys.argv[1:]:
  if not arg.startswith('--'): continue
  toks = arg.split('=')
  key = toks.pop(0)[2:]
  value = (not len(toks)) and True or '='.join(toks)
  module_kwargs[key] = value

if "__main__" == __name__:
  tt = TILTED_TANK(**module_kwargs)
  print((tt.volcalc(tt.upperlim), tt.volest(tt.upperlim),))
  print(tt)

def calc_vs_est(R=1.0, phi=30.0, **kwargs):

  import matplotlib.pyplot as plt

  log2s = range(7,15,1)

  ns = [1<<log2 for log2 in range(6,16)]
  tts = [TILTED_TANK(R=R, phi=phi, n=n) for n in ns]

  Level100 = tts[0].upperlim / 100.0

  percents = [percent > 0 and percent or 1 for percent in range(0,101,10)]

  Volumes = list()

  for percent in percents:
    Level = percent * Level100
    volcalc = tts[0].volcalc(Level)
    Volumes.append(volcalc)
    fracdiffs = [(abs(tt.volest(Level) - volcalc) / volcalc) for tt in tts]
    plt.loglog(ns, fracdiffs, label=f"{percent}%")

  plt.title(f"Numerical vs. Analytical Volume Models\n{tts[-1]}")
  plt.xlabel('Numerical model cell count across diameter, n')
  plt.ylabel('Fractional volume difference, Numerical vs. Analytical models')
  plt.legend()
  plt.show()

  plt.plot(percents,Volumes)
  plt.title(f"Analytical Volume Model\n{tts[-1]} (ignore n)")
  plt.xlabel('Level, %')
  plt.ylabel('Volume, $R^{3}$ units')
  plt.show()
  
