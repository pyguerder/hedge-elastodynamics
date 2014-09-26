# -*- coding: utf-8 -*-
"""Functions used for generating the sources"""

from __future__ import division

__copyright__ = "Copyright (C) 2011 Olivier Bou Matar <olivier.boumatar@iemn.univ-lille1.fr>"

__license__ = """
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see U{http://www.gnu.org/licenses/}.
"""

from hedge.data import ITimeDependentGivenFunction
from math import sin, pi, exp


class TimeRickerWaveletGivenFunction(ITimeDependentGivenFunction):
    """Modulates an :class:`ITimeDependentGivenFunction` by a Ricker Wavelet
    in time.
    """
    def __init__(self, gf, fc, tD):
        self.gf = gf
        self.fc = fc
        self.tD = tD

    def volume_interpolant(self, t, discr):
        return (0.5 - (pi*self.fc)**2 * (t - self.tD)**2) * exp(-(pi*self.fc)**2 * (t - self.tD)**2)\
                * self.gf.volume_interpolant(t, discr)

    def boundary_interpolant(self, t, discr, tag):
        return (0.5 - (pi*self.fc)**2 * (t - self.tD)**2) * exp(-(pi*self.fc)**2 * (t - self.tD)**2)\
                * self.gf.boundary_interpolant(t, discr, tag)


class SinusGivenFunction(ITimeDependentGivenFunction):
    """Modulates an :class:`ITimeDependentGivenFunction` by a sinus in time.
    """
    def __init__(self, gf, fc, tD):
        self.gf = gf
        self.fc = fc
        self.tD = tD

    def volume_interpolant(self, t, discr):
        return (0.5 * sin(2*pi*self.fc*t))\
                * self.gf.volume_interpolant(t, discr)

    def boundary_interpolant(self, t, discr, tag):
        return (0.5 * sin(2*pi*self.fc*t))\
                * self.gf.boundary_interpolant(t, discr, tag)

class ModulatedSinusGivenFunction(ITimeDependentGivenFunction):
    """Modulates an :class:`ITimeDependentGivenFunction` by an amplitude modulated sinus in time.
    """
    def __init__(self, gf, fc, im):
        self.gf = gf
        self.fc = fc
        self.im = im

    def volume_interpolant(self, t, discr):
        return 0.5 * (sin(2*pi*self.fc*t) * self.im + 1.0)\
                * self.gf.volume_interpolant(t, discr)

    def boundary_interpolant(self, t, discr, tag):
        return 0.5 * (sin(2*pi*self.fc*t) * self.im + 1.0)\
                * self.gf.boundary_interpolant(t, discr, tag)

class SineBurstGivenFunction(ITimeDependentGivenFunction):
    """Modulates an :class:`ITimeDependentGivenFunction` by a sinus in an envelope in time.
    """
    def __init__(self, gf, fc, tD):
        self.gf = gf
        self.fc = fc
        self.tD = tD

    def volume_interpolant(self, t, discr):
        return 5e4 * sin(2 * pi * self.fc * t) * exp( - ((t - 10 / self.fc) * self.fc / 4)**2) \
                * self.gf.volume_interpolant(t, discr)

    def boundary_interpolant(self, t, discr, tag):
        return 5e4 * sin(2 * pi * self.fc * t) * exp( - ((t - 10 / self.fc) * self.fc / 4)**2) \
                * self.gf.boundary_interpolant(t, discr, tag)
