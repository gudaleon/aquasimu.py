"""Port target of aqcore/sedcomp.cpp (SEDCOMP::Delta, GridValue, SetGridPts, FluxLimit).

1-D sediment column, z from Zt (surface) to Zb, numgrid points (equal cells, dz=(Zb-Zt)/(n-2)).
Per cell unknowns: Q (volumetric solid flux), porosity theta, then state vars.
Dissolved C:  d(theta C)/dt = -d(q_w C)/dz + d/dz(theta D_s dC/dz) + theta r_C
Particulate X: d((1-theta) X)/dt = -d(q_s X)/dz + d/dz(D_b d X/dz) + (1-theta) r_X
Porosity: d theta/dt = RatePor (user formula), surface porosity SurfPor.
Top BC: flux (particulates) or concentration (dissolved); bottom BC: concentration (2003 revision).
Advection uses flux limiter (FluxLimit) when HighRes=TRUE, else upwind.
TODO v1: implement residual() exactly per sedcomp.cpp lines ~1203-1845.
"""
import numpy as np


class SedComp:
    def __init__(self, numgrid=250, zt=0.0, zb=0.19, highres=False):
        self.n, self.zt, self.zb, self.highres = numgrid, zt, zb, highres

    @property
    def dz(self):
        return (self.zb - self.zt) / (self.n - 2)

    def grid(self):
        z = np.empty(self.n)
        z[0] = self.zt; z[-1] = self.zb
        z[1:-1] = self.zt + self.dz * (np.arange(self.n - 2) + 0.5)
        return z

    def residual(self, t, y, ydot, params):
        raise NotImplementedError("port SEDCOMP::Delta")
