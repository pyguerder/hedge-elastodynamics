# -*- coding: utf-8 -*-
"""Functions for converting the matrices and tensors and converting the indices"""

from __future__ import division

__copyright__ = "Copyright (C) 2011 Pierre-Yves Guerder <pierre-yves.guerder@centraliens-lille.org>"

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

import numpy

class Utils:

    conv_2D = numpy.array([
                    [1, 3],
                    [4, 2]])
    conv_3D = numpy.array([
                    [1, 6, 5],
                    [9, 2, 4],
                    [8, 7, 3]
                    ])
    conv_2D_sym = numpy.array([
                    [1, 6],
                    [6, 2]])
    conv_3D_sym = numpy.array([
                    [1, 6, 5],
                    [6, 2, 4],
                    [5, 4, 3]
                    ])


    def convert_dim(C, dim):

        conv = numpy.zeros(dim**2, dtype=int)

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
            raise RuntimeError, "Bad number of dimensions"

        Ceq = numpy.zeros((dim**2, dim**2), dtype=C.dtype)
        for i in range(dim**2):
            for j in range(dim**2):
                Ceq[i,j] = C[conv[i]-1,conv[j]-1]
        return Ceq

    def condense(i, j, dim):
        if dim == 1:
            return 1
        elif dim == 2:
            return Utils.conv_2D[i, j]
        elif dim == 3:
            return Utils.conv_3D[i, j]

    def condense_sym(i, j, dim):
        if dim == 1:
            return 1
        elif dim == 2:
            return Utils.conv_2D_sym[i, j]
        elif dim == 3:
            return Utils.conv_3D_sym[i, j]

    def kronecker(i, j):
        if i == j :
            return 1
        else :
            return 0

    convert_dim = staticmethod(convert_dim)
    kronecker = staticmethod(kronecker)
    condense = staticmethod(condense)
    condense_sym = staticmethod(condense_sym)
