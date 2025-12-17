"""
Model volume of liquid in a flat-bottomed cylindrical tank resting on a
sloped base so the tank's cylindrical axis is tilted at angle phi from
vertical and the flat bottom is at angle phi from horizontal.

Usage:

  python tilted_tank.py [--R=<radius>] [--phi=<tilt,deg>] [--n=<row count>]

"""
import os
import sys
import numpy

########################################################################
class TILTED_TANK:
  """
Cylindrical tank model class
* Radius R, diameter 2R
* Flat bottom is tilted at angle phi wrt horizondal
* Models between upper and lower level limits:
* 0 = Lower limit => low side of bottom
* 2R sin(phi) = Upper limit => high side of bottom


                              /
                             /
                            /
                           /     Liquid
                    |     /      surface      /
                    |    /          |        /  |
Highest point of    v   /           |       /   |
bottom of tank => ----- --__        v      /    v
                       |    --__==========/ -------
Tilt angle phi =>      |        --__     /             
                  _____|__________  --__/  ________ <= Lowest level
                    ^                           ^       in tank,
                    |                           |       value = 0
                    |                           |
     Highest measured                           Level measurement
      level of model,                           distance above lowest
  value = 2R sin(phi)                           point in tank

"""


  ######################################################################
  def __init__(self, R=1.0, phideg=30, n=64, **kwargs):
    """
  R - radius of tank
phi - angle of tilt, degrees
  n - number of rows across the diameter (2R) of the tank

Numerically model the volume of liquid partially covering the tilted
flat bottom of a cylindrical tank using [n] rows with their long axis
parallel to the rotation axis of the tilt, and with all rows arranged
side-by-side in a plane parallel to the tilted bottom of the tank,
where each row
- has a square cross-section, [2R/n wide] x [2R/n high], and
- has a length of 2sqrt(R^2 - P^2), so it fits inside the tank walls @ P
where P is the offset of the centerline of the row from the tank's axis

"""
    ### Add input parameters to model;
    ### Convert angle argument phi units to radians in self.phi
    self.R, self.phideg, self.n = float(R), float(phideg), int(n)
    self.phi = self.phideg * numpy.pi / 180.0

    ### Calculate

    ### - trigonometric values from angle phi
    self.cos,self.sin = numpy.cos(self.phi), numpy.sin(self.phi)
    self.tan = numpy.tan(self.phi)

    ### - model upper limit for level i.e. high side of tank bottom
    self.upperlim = 2.0 * self.R * self.sin

    ### - row width- and height-related values
    self.halfedge = self.R / self.n       ### Offset to row center
    self.edge = 2.0 * self.halfedge       ### Row width and height

    ### - Offsets of rows' centerlines from low point of tank bottom
    self.rowxs = (self.edge * numpy.arange(self.n)) + self.halfedge

    ### - Lengths of rows, Y = sqrt(R^2 - P^2), where P = R-x
    ### - Volumes of rows
    ylens = 2.0 * numpy.sqrt((self.R**2) - ((self.R-self.rowxs)**2))
    self.rowvols = (self.edge**2) * ylens

    ### - Vertical offset of tilted (rotated) row bottoms wrt level=0
    self.rowbaseys = self.rowxs * self.sin

    ### - Vertical offset between successive tilted row bottoms
    self.rowdy = self.edge * self.cos


  ######################################################################
  def __repr__(self):
    return f"<{self.__class__.__name__} R:{self.R} phi:{self.phideg}deg n:{self.n}>"


  ######################################################################
  def volest(self, Level):
    """
Numerically estimate the volume of liquid partially covering the tilted
flat bottom of a cylindrical tank

Instead of incrementing the elevation of the rows, decrement the level
against which the rows' levels are compared

"""

    ### Clamp the level argument to the valid range, and reduce it
    ### by one-half of the row-to-row offset
    Llim = max([0.0,min([Level,self.upperlim])]) - (self.rowdy / 2.0)

    ### Initialize the volume estimate
    Vest = 0.0

    ### Loop while the current level value is positive
    while Llim > 0.0:
      ### Select rows that are below the current level
      ### If no rows are selected, exit the loop
      ### Sum the volumes of the selected rows to the volume estimate
      ### Decrement the level value by one rotated row height
      iw = numpy.where(self.rowbaseys < Llim)
      if not len(iw[0]): break
      Vest += numpy.sum(self.rowvols[iw])
      Llim -= self.rowdy

    return Vest


  ######################################################################
  def volcalc(self, Level):
    """
Analytically calculate the volume of liquid partially covering the tilted
flat bottom of a cylindrical tank

"""
    ### H - height of wet segment on tilted circular flat bottom
    ### X - cosine of half-angle of wet segment

    H = Level / self.sin
    X = (self.R - H) / self.R

    ### Limit X value if necessary, and adjust H
    if X   >  1.0: X,H =  1.0,0.0
    elif X < -1.0: X,H = -1.0,2.0 * self.R

    ### Calculate (2RH - H^2) parameter used in segment area formula
    Y = (2.0 * self.R * H) - (H**2)

    ### Calculate volume
    Xterm = (self.R**3) * (numpy.sqrt(1.0 - (X**2)) - (X*numpy.arccos(X)))
    Yterm = - (Y**1.5) / 3.0
    return self.tan * (Xterm + Yterm)


########################################################################
module_kwargs = dict()
for arg in sys.argv[1:]:
  if not arg.startswith('--'): continue
  toks = arg.split('=')
  key = toks.pop(0)[2:]
  value = (not len(toks)) and True or '='.join(toks)
  module_kwargs[key] = value


########################################################################
if "__main__" == __name__:
  tt = TILTED_TANK(**module_kwargs)
  print((tt.volcalc(tt.upperlim), tt.volest(tt.upperlim),))
  print(tt)


########################################################################
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
  
