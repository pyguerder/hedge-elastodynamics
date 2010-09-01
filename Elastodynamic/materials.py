# Hedge - the Hybrid'n'Easy DG Environment
# Copyright (C) 2007 Andreas Kloeckner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import numpy

def Properties_Mat(dtype,material="Aluminium",dim=2):
    
    CIJ = numpy.zeros((dim**2,dim**2),dtype=dtype)
    C = numpy.zeros((7,7),dtype=dtype)
    conv = numpy.zeros(dim**2,dtype=int)

    if dim == 1:
        conv[0] = 1
    elif dim == 2:
        conv[0] = 1
        conv[1] = 2
        conv[2] = 6
        conv[3] = 6
    elif dim == 3:
        conv[0] = 1
        conv[1] = 2
        conv[2] = 3
        conv[3] = 4
        conv[4] = 5
        conv[5] = 6
        conv[6] = 4
        conv[7] = 5
        conv[8] = 6
    else:
        raise RuntimeError, "bad number of dimensions"
    
    if material == "Aluminium":
	rho = 2700.
        C[1,1] = rho*6410.0**2
        C[1,2] = rho*(6410.0**2-2*2907.0**2)
        C[1,3] = C[1,2]
        C[1,4] = 0.0
        C[1,5] = 0.0
        C[1,6] = 0.0

        C[2,1] = C[1,2]
        C[2,2] = C[1,1]
        C[2,3] = C[1,2]
        C[2,4] = 0.0
        C[2,5] = 0.0
        C[2,6] = 0.0

        C[3,1] = C[1,2]
        C[3,2] = C[1,2]
        C[3,3] = C[1,1]
        C[3,4] = 0.0
        C[3,5] = 0.0
        C[3,6] = 0.0

        C[4,1] = 0.0
        C[4,2] = 0.0
        C[4,3] = 0.0
        C[4,4] = rho*2907.0**2
        C[4,5] = 0.0
        C[4,6] = 0.0

        C[5,1] = 0.0
        C[5,2] = 0.0
        C[5,3] = 0.0
        C[5,4] = 0.0
        C[5,5] = C[4,4]
        C[5,6] = 0.0

        C[6,1] = 0.0
        C[6,2] = 0.0
        C[6,3] = 0.0
        C[6,4] = 0.0
        C[6,5] = 0.0
        C[6,6] = C[4,4]
    else:
	raise RuntimeError, "Material not defined"

    for i in range(dim**2):
	for j in range(dim**2):
            CIJ[i,j] = C[conv[i],conv[j]]    

    return rho, CIJ
