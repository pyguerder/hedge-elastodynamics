# -*- coding: utf8 -*-
"""Operator for non linear Elastodynamics equations."""

from __future__ import division

__copyright__ = "Copyright (C) 2007 Hendrik Riedmann, Andreas Kloeckner"

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
import hedge.tools
import hedge.mesh
import hedge.data
from hedge.models import HyperbolicOperator
from pytools import Record
from hedge.tools import is_zero

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

        Ceq = numpy.zeros((dim**2, dim**2))
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


class NLElastoDynamicsOperator(HyperbolicOperator):
    """An nD non linear Elastodynamics operator.

    see YiFeng LI PhD p. 41

    dq/dt = -dF/dx  - dG/dy - dH/dz

    where e.g. in 3D

    q = (rho_v_1, rho_v_2, rho_v_3, F_11, F_22, F_33, F_23, F_13, F_12, F_32, F_31, F_21)
    F = (-P_11, -P_21, -P_31, -v_1, 0, 0, 0, 0, 0, 0, -v_3, -v_2)
    G = (-P_12, -P_22, -P_32, 0, -v_2, 0, 0, 0, -v_1, -v_3, 0, 0)
    H = (-P_13, -P_23, -P_33, 0, 0, -v_3, -v_2, -v_1, 0, 0, 0, 0)

    For the linear case with included attenuation
    the stress-strain relation is given by

    P_ij = C_ijkl * F_kl + nu_ijkl * dv_k/dx_l

    Field order is [rho_v_1, rho_v_2, rho_v_3, F_11, F_22, F_33, F_23, F_13, F_12, F_32, F_31, F_21].
    """
    def __init__(self,
                 dimensions,
                 rho,
                 C,                             # 6×6 matrix
                 Cnl,                           # 6×6×6 tensor
                 nonlinearity_type = "classical",
                 stressfree_tag="stressfree",
                 fixed_tag="fixed",
                 open_tag="open",
                 source=None,
                 flux_type = "lf",
            ):
        """
        :param source: should implement
        :class:`hedge.data.IFieldDependentGivenFunction`
        or be None.
        """

        self.dimensions = dimensions

        self.C = C
        self.Ceq = Utils.convert_dim(self.C, self.dimensions)
        self.Cnl = Cnl

        self.MIJK = self.M()
        self.rho = rho

        self.stressfree_tag = stressfree_tag
        self.fixed_tag = fixed_tag
        self.open_tag = open_tag

        self.source = source
        self.flux_type = flux_type

        self.nonlinearity_type = nonlinearity_type

    def rho_v(self, state):
        return state[0:self.dimensions]

    def F(self, state):
        return state[self.dimensions:(self.dimensions+1)*self.dimensions]

    def v(self, state):
        from pytools.obj_array import make_obj_array
        return make_obj_array([
                rho_v_i/self.rho
                for rho_v_i in self.rho_v(state)])

    def P(self, state):
        Pi = numpy.zeros((self.dimensions**2), dtype=object)
        Ceq = self.C_eq(state)
        F = self.F(state)
        for i in range(self.dimensions**2):
            for j in range(self.dimensions**2):
                Pi[i] += Ceq[i,j]*F[j]
        return Pi

    def M(self):
        dim = self.dimensions
        M = numpy.zeros((dim**2, dim**2, dim**2), dtype=object)
        for i in range(dim):
            for j in range(dim):
                I1 = Utils.condense(i, j, dim) - 1
                I2 = Utils.condense_sym(i, j, dim) - 1
                for k in range(dim):
                    for l in range(dim):
                        J1 = Utils.condense(k, l, dim) - 1
                        J2 = Utils.condense_sym(k, l, dim) - 1
                        for m in range(dim):
                            for n in range(dim):
                                K1 = Utils.condense(m, n, dim) - 1
                                K2 = Utils.condense_sym(m, n, dim) - 1
                                LM = Utils.condense_sym(l, m, dim) - 1
                                IL = Utils.condense_sym(i, l, dim) - 1
                                IK = Utils.condense_sym(i, k, dim) - 1
                                LM = Utils.condense_sym(l, m, dim) - 1
                                M[I1,J1,K1] = self.Cnl[I2,J2,K2] \
                                    + self.C[I2, LM] * Utils.kronecker(k,n) \
                                    + self.C[IL, K2] * Utils.kronecker(j,k) \
                                    + self.C[IK, LM] * Utils.kronecker(j,n)
        return M

    def C_eq(self, state):
        if self.nonlinearity_type == "classical":
            C = numpy.zeros((self.dimensions**2,self.dimensions**2), dtype=object)
            F = self.F(state)
            for i in range(self.dimensions**2):
                for j in range(self.dimensions**2):
                    for k in range(self.dimensions**2):
                        C[i,j] += self.Ceq[i,j] + 0.5 * self.MIJK[i, j, k] * F[k]
        else:
            C = self.Ceq
        return C

    def flux_num(self, state, fluxes, bdry_tags_state_and_fluxes):
        from hedge.flux import FluxVectorPlaceholder, make_normal
        from hedge.optemplate import get_flux_operator, BoundaryPair
        from pytools.obj_array import join_fields
        from math import sqrt

        n = len(state)
        d = len(fluxes)
        fvph = FluxVectorPlaceholder(n*(d+1))
        state_ph  = fvph[0:n]
        fluxes_ph = [fvph[i*n:(i+1)*n] for i in range(1,d+1)]

        normal = make_normal(d)
        penalty = sqrt(self.C[0,0]/self.rho)*(state_ph.ext-state_ph.int)

        flux_strong = 0.5*sum(n_i*(f_i.int-f_i.ext) for n_i, f_i in zip(normal, fluxes_ph))

        if self.flux_type == "central":
            pass
        elif self.flux_type == "lf":
            flux_strong = 0.5*penalty + flux_strong
        else:
            raise ValueError("invalid flux type '%s'" % self.flux_type)

        flux_op = get_flux_operator(flux_strong)
        int_operand = join_fields(state,*fluxes)

        return (flux_op(int_operand)
                +sum(flux_op(BoundaryPair(int_operand, join_fields(bdry_state, *bdry_fluxes), tag))
                        for tag, bdry_state, bdry_fluxes in bdry_tags_state_and_fluxes))

    def op_template(self):
        from hedge.optemplate import \
                Field, \
                BoundaryPair, \
                get_flux_operator, \
                make_nabla, \
                InverseMassOperator, \
                BoundarizeOperator
        from hedge.optemplate.tools import make_vector_field

        from pytools.obj_array import join_fields
        from hedge.tools.symbolic import make_common_subexpression as cse

        dim = self.dimensions

        def flux(q):
            if dim == 1:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[0],

                            # flux F
                            -self.v(q)[0]
                            ), "x_flux")
                        ]

            elif dim == 2:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[0],-self.P(q)[3],

                            # flux F
                            -self.v(q)[0],0,0,-self.v(q)[1]
                            ), "x_flux"),
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[2],-self.P(q)[1],

                            # flux F
                            0,-self.v(q)[1],-self.v(q)[0],0
                            ), "y_flux")
                        ]
            else:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[0],-self.P(q)[8],-self.P(q)[7],

                            # flux F
                            -self.v(q)[0],0,0,0,0,0,0,-self.v(q)[2],-self.v(q)[1]
                            ), "x_flux"),
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[5],-self.P(q)[1],-self.P(q)[6],

                            # flux F
                            0,-self.v(q)[1],0,0,0,-self.v(q)[0],-self.v(q)[2],0,0
                            ), "y_flux"),
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[4],-self.P(q)[3],-self.P(q)[2],

                            # flux F
                            0,0,-self.v(q)[2],-self.v(q)[1],-self.v(q)[0],0,0,0,0
                            ), "z_flux")
                        ]

        def bdry_flux(q_bdry, tag):
            # stress free BCs -------------------------------------------------------
            if tag == self.stressfree_tag:
                signP = 1
                signv = -1
            # fixed BCs -------------------------------------------------------------
            elif tag == self.fixed_tag:
                signP = -1
                signv = 1
            else:
                raise ValueError("invalid Bdry conditions")

            if dim == 1:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[0],

                            # flux F
                            signv*self.v(q_bdry)[0]
                            ), "x_bflux")
                        ]
            elif dim == 2:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[0],signP*self.P(q_bdry)[3],

                            # flux F
                            signv*self.v(q_bdry)[0],0,0,signv*self.v(q_bdry)[1]
                            ), "x_bflux"),
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[2],signP*self.P(q_bdry)[1],

                            # flux F
                            0,signv*self.v(q_bdry)[1],signv*self.v(q_bdry)[0],0
                            ), "y_bflux")
                        ]
            else:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[0],signP*self.P(q_bdry)[8],signP*self.P(q_bdry)[7],

                            # flux F
                            signv*self.v(q_bdry)[0],0,0,0,0,0,0,signv*self.v(q_bdry)[2],signv*self.v(q_bdry)[1]
                            ), "x_bflux"),
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[5],signP*self.P(q_bdry)[1],signP*self.P(q_bdry)[6],

                            # flux F
                            0,signv*self.v(q_bdry)[1],0,0,0,signv*self.v(q_bdry)[0],signv*self.v(q_bdry)[2],0,0
                            ), "y_bflux"),
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[4],signP*self.P(q_bdry)[3],signP*self.P(q_bdry)[2],

                            # flux F
                            0,0,signv*self.v(q_bdry)[2],signv*self.v(q_bdry)[1],signv*self.v(q_bdry)[0],0,0,0,0
                            ), "z_bflux"),
                        ]


        state = make_vector_field("q", dim*(dim+1))
        #rho_v = state[0:dim]
        #F = state[dim:dim*(dim+1)]

        # fluxes ------------------------------------------------------------
        fluxes = flux(state)

        # boundary conditions ---------------------------------------------------
        state_bc_stressfree = BoundarizeOperator(self.stressfree_tag)(state)
        state_bc_fixed      = BoundarizeOperator(self.fixed_tag)(state)

        all_tags_and_bcs = [
                (self.stressfree_tag, state_bc_stressfree),
                (self.fixed_tag, state_bc_fixed)
                           ]

        bdry_tags_state_and_fluxes = [(tag, bc, bdry_flux(bc, tag)) for tag, bc in all_tags_and_bcs]

        # entire operator -----------------------------------------------------
        from math import sin, cos, pi
        nabla = make_nabla(dim)

        result = (-numpy.dot(nabla,fluxes) + InverseMassOperator() * (self.flux_num(state,fluxes,bdry_tags_state_and_fluxes)))

        if self.source is not None:
            result[0] += Field("source_v_x")*sin(10*pi/180)
            result[1] += Field("source_v_x")*cos(10*pi/180)

        return result

    def bind(self, discr):
        from hedge.mesh import check_bc_coverage
        check_bc_coverage(discr.mesh, [
            self.stressfree_tag,
            self.fixed_tag,
            self.open_tag])

        compiled_op_template = discr.compile(self.op_template())

        def rhs(t, q):
            extra_kwargs = {}

            if self.source is not None:
                extra_kwargs["source_v_x"] = self.source.volume_interpolant(t, discr)

            return compiled_op_template(q=q,**extra_kwargs)

        return rhs

    def max_eigenvalue(self, t, fields=None, discr=None):
        from math import sqrt

        return sqrt(self.C[0,0]/self.rho)



class ElastoDynamicsOperator(HyperbolicOperator):
    """An nD linear Elastodynamics operator.

    dq/dt = -dF/dx  - dG/dy - dH/dz

    where e.g. in 3D

    q = (rho_v_1, rho_v_2, rho_v_3, F_11, F_22, F_33, 2*F_23, 2*F_13, 2*F_12)
    F = (-P_11, -P_12, -P_13, -v_1, 0, 0, 0, -v_3, -v_2)
    G = (-P_12, -P_22, -P_23, 0, -v_2, 0, -v_3, 0, -v_1)
    H = (-P_13, -P_23, -P_33, 0, 0, -v_3, -v_2, -v_1, 0)

    For the linear case with included attenuation
    the stress-strain relation is given by

    P_ij = C_ijkl * F_kl + nu_ijkl * dv_k/dx_l

    Field order is [rho_v_1, rho_v_2, rho_v_3, F_11, F_22, F_33, 2*F_23, 2*F_13, 2*F_12].
    """
    def __init__(self, dimensions, rho,
            C, stressfree_tag="stressfree",
            fixed_tag="fixed",
            open_tag="open",
            source=None,
            flux_type = "lf",
            dimF=None,
            ):
        """
        :param source: should implement
        :class:`hedge.data.IFieldDependentGivenFunction`
        or be None.
        """

        self.dimensions = dimensions
        self.dimF = [0, 1, 3, 6]

        self.C = C
        self.Ceq = Utils.convert_dim(self.C, self.dimensions)
        self.rho = rho

        self.stressfree_tag = stressfree_tag
        self.fixed_tag = fixed_tag
        self.open_tag = open_tag

        self.source = source
        self.flux_type = flux_type

    def rho_v(self, q):
        return q[0:self.dimensions]

    def F(self, q):
        return q[self.dimensions:self.dimensions+self.dimF[self.dimensions]]

    def v(self, q):
        from pytools.obj_array import make_obj_array
        return make_obj_array([
                rho_v_i/self.rho
                for rho_v_i in self.rho_v(q)])

    def P(self, q):
        Pi = numpy.zeros((self.dimF[self.dimensions]), dtype=object)
        F = self.F(q)
        for i in range(self.dimF[self.dimensions]):
            for j in range(self.dimF[self.dimensions]):
                    Pi[i] += self.Ceq[i,j]*F[j]
        return Pi

    def cse_rho_v(self, q):
        from hedge.tools.symbolic import make_common_subexpression as cse
        return cse(self.rho_v(q), "rho_v")

    def cse_F(self, q):
        from hedge.tools.symbolic import make_common_subexpression as cse
        return cse(self.F(q), "F")

    def cse_v(self, q):
        from hedge.tools.symbolic import make_common_subexpression as cse
        return cse(self.v(q), "v")

    def cse_P(self, q):
        from hedge.tools.symbolic import make_common_subexpression as cse
        return cse(self.P(q), "P")

    def flux_num(self, q, fluxes, bdry_tags_state_and_fluxes):
        from hedge.flux import FluxVectorPlaceholder, make_normal
        from hedge.optemplate import get_flux_operator, BoundaryPair

        from pytools.obj_array import join_fields
        from math import sqrt

        n = len(q)
        d = len(fluxes)
        fvph = FluxVectorPlaceholder(n*(d+1))
        state_ph  = fvph[0:n]
        fluxes_ph = [fvph[i*n:(i+1)*n] for i in range(1,d+1)]

        normal = make_normal(d)
        penalty = sqrt(self.Ceq[0,0]/self.rho)*(state_ph.ext-state_ph.int)

        flux_strong = 0.5*sum(n_i*(f_i.int-f_i.ext) for n_i, f_i in zip(normal, fluxes_ph))

        if self.flux_type == "central":
            pass
        elif self.flux_type == "lf":

            flux_strong = 0.5*penalty + flux_strong
        else:
            raise ValueError("invalid flux type '%s'" % self.flux_type)

        flux_op = get_flux_operator(flux_strong)
        int_operand = join_fields(q,*fluxes)

        return (flux_op(int_operand)
                +sum(flux_op(BoundaryPair(int_operand, join_fields(bdry_state, *bdry_fluxes), tag))
                      for tag, bdry_state, bdry_fluxes in bdry_tags_state_and_fluxes))

    def op_template(self):
        from hedge.optemplate import \
                Field, \
                BoundaryPair, \
                get_flux_operator, \
                make_nabla, \
                InverseMassOperator, \
                BoundarizeOperator
        from hedge.optemplate.tools import make_vector_field

        from pytools.obj_array import join_fields
        from hedge.tools.symbolic import make_common_subexpression as cse

        dim = self.dimensions

        def flux(q):
            if dim == 1:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[0],

                            # flux F
                            -self.v(q)[0]
                            ), "x_flux")
                        ]

            elif dim == 2:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[0],-self.P(q)[2],

                            # flux F
                            -self.v(q)[0],0,-self.v(q)[1]
                            ), "x_flux"),
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[2],-self.P(q)[1],

                            # flux F
                            0,-self.v(q)[1],-self.v(q)[0]
                            ), "y_flux")
                        ]
            elif dim == 3:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[0],-self.P(q)[5],-self.P(q)[4],

                            # flux F
                            -self.v(q)[0],0,0,0,-self.v(q)[2],-self.v(q)[1]
                            ), "x_flux"),
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[5],-self.P(q)[1],-self.P(q)[3],

                            # flux F
                            0,-self.v(q)[1],0,-self.v(q)[2],0,-self.v(q)[0]
                            ), "y_flux"),
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[4],-self.P(q)[3],-self.P(q)[2],

                            # flux F
                            0,0,-self.v(q)[2],-self.v(q)[1],-self.v(q)[0],0
                            ), "z_flux")
                        ]
            else:
                raise ValueError("invalid dimension")

        def bdry_flux(q_bdry, tag):
            # stress free BCs -------------------------------------------------------
            if tag == self.stressfree_tag:
                signP = 1
                signv = -1
            # fixed BCs -------------------------------------------------------------
            elif tag == self.fixed_tag:
                signP = -1
                signv = 1
            else:
                raise ValueError("invalid Bdry conditions")

            if dim == 1:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[0],

                            # flux F
                            signv*self.v(q_bdry)[0]
                            ), "x_bflux")
                        ]
            elif dim == 2:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[0],signP*self.P(q_bdry)[2],

                            # flux F
                            signv*self.v(q_bdry)[0],0,signv*self.v(q_bdry)[1]
                            ), "x_bflux"),
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[2],signP*self.P(q_bdry)[1],

                            # flux F
                            0,signv*self.v(q_bdry)[1],signv*self.v(q_bdry)[0]
                            ), "y_bflux")
                        ]
            elif dim == 3:
                return [ # one entry for each flux direction
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[0],signP*self.P(q_bdry)[5],signP*self.P(q_bdry)[4],

                            # flux F
                            signv*self.v(q_bdry)[0],0,0,0,signv*self.v(q_bdry)[2],signv*self.v(q_bdry)[1]
                            ), "x_bflux"),
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[5],signP*self.P(q_bdry)[1],signP*self.P(q_bdry)[3],

                            # flux F
                            0,signv*self.v(q_bdry)[1],0,signv*self.v(q_bdry)[2],0,signv*self.v(q_bdry)[0]
                            ), "y_bflux"),
                        cse(join_fields(
                            # flux rho_v
                            signP*self.P(q_bdry)[4],signP*self.P(q_bdry)[3],signP*self.P(q_bdry)[2],

                            # flux F
                            0,0,signv*self.v(q_bdry)[2],signv*self.v(q_bdry)[1],signv*self.v(q_bdry)[0],0
                            ), "z_bflux")
                        ]
            else:
                raise ValueError("invalid dimension")


        state = make_vector_field("q", dim+self.dimF[dim])

        # fluxes ------------------------------------------------------------
        fluxes = flux(state)

        # boundary conditions ---------------------------------------------------
        state_bc_stressfree = BoundarizeOperator(self.stressfree_tag)(state)
        state_bc_fixed      = BoundarizeOperator(self.fixed_tag)(state)

        all_tags_and_bcs = [
                (self.stressfree_tag, state_bc_stressfree),
                (self.fixed_tag, state_bc_fixed)
                           ]

        bdry_tags_state_and_fluxes = [(tag, bc, bdry_flux(bc, tag)) for tag, bc in all_tags_and_bcs]

        # entire operator -----------------------------------------------------
        from math import sin, cos, pi
        nabla = make_nabla(dim)

        result = (-numpy.dot(nabla,fluxes) + InverseMassOperator() * (self.flux_num(state,fluxes,bdry_tags_state_and_fluxes)))

        if self.source is not None:
            result[0] += Field("source_v_x")*sin(10*pi/180)
            result[1] += Field("source_v_x")*cos(10*pi/180)

        return result

    def bind(self, discr):
        from hedge.mesh import check_bc_coverage
        check_bc_coverage(discr.mesh, [
            self.stressfree_tag,
            self.fixed_tag,
            self.open_tag])

        compiled_op_template = discr.compile(self.op_template())

        def rhs(t, q):
            extra_kwargs = {}

            if self.source is not None:
                extra_kwargs["source_v_x"] = self.source.volume_interpolant(t, discr)

            return compiled_op_template(q=q,**extra_kwargs)

        return rhs

    def max_eigenvalue(self, t, fields=None, discr=None):
        from math import sqrt

        return sqrt(self.Ceq[0,0]/self.rho)




# To be finished !!!!
class NPMLElastoDynamicsOperator(ElastoDynamicsOperator):
    """Implements a NPML as in

    YiFeng LI Ph'D p. 121-138
    """

    class PMLCoefficients(Record):
        __slots__ = ["sigma", "alpha", "kappa"]

        def map(self, f):
            return self.__class__(
                    **dict((name, f(getattr(self, name)))
                        for name in self.fields))

    def __init__(self, *args, **kwargs):
        #self.add_decay = kwargs.pop("add_decay", True)
        ElastoDynamicsOperator.__init__(self, *args, **kwargs)

    def pml_local_op(self, w):
        dim = self.dimensions

        from hedge.optemplate import make_vector_field
        sigma = make_vector_field("sigma", self.dimensions)
        alpha = make_vector_field("alpha", self.dimensions)
        kappa = make_vector_field("kappa", self.dimensions)


        rhs = numpy.zeros((2*self.dimensions**2,), dtype=object)

        #for mx in range(3):
        #    my = (mx+1) % 3
        #    mz = (mx+2) % 3

        #    from hedge.tools.mathematics import levi_civita
        #    assert levi_civita((mx,my,mz)) == 1

        #    rhs[mx] += -sig[my]/self.epsilon*(2*e[mx]+p[mx]) - 2*tau[my]/self.epsilon*e[mx]
        #    rhs[my] += -sig[mx]/self.epsilon*(2*e[my]+p[my]) - 2*tau[mx]/self.epsilon*e[my]
        #    rhs[3+mz] += 1/(self.epsilon*self.mu) * (
        #      sig_prime[mx] * q[mx] - sig_prime[my] * q[my])

        #    rhs[6+mx] += sig[my]/self.epsilon*e[mx]
        #    rhs[6+my] += sig[mx]/self.epsilon*e[my]
        #    rhs[9+mx] += -sig[mx]/self.epsilon*q[mx] - (e[my] + e[mz])

        return rhs

    def op_template(self, q=None):
        if q is None:
            from hedge.optemplate import make_vector_field
            q = make_vector_field("q", dim+self.dimF[dim])

        from hedge.tools import join_fields
        return join_fields(
                ElastoDynamicsOperator.op_template(self, q),
                numpy.zeros((2*self.dimensions**2,), dtype=object)
                ) + self.pml_local_op(q)

    def bind(self, discr, coefficients):
        return ElastoDynamicsOperator.bind(self, discr,
                sigma=coefficients.sigma,
                alpha=coefficients.alpha,
                kappa=coefficients.kappa)


    # sigma business ----------------------------------------------------------
    def _construct_scalar_coefficients(self, discr, node_coord,
            i_min, i_max, o_min, o_max, sigma_exponent, alpha_exponent, kappa_exponent):
        assert o_min < i_min <= i_max < o_max

        if o_min != i_min:
            l_dist = (i_min - node_coord) / (i_min-o_min)
            l_dist[l_dist < 0] = 0.0
            lp_dist = 1.0-(i_min - node_coord) / (i_min-o_min)
            lp_dist[l_dist < 0] = 0.0
        else:
            l_dist = lp_dist = numpy.zeros_like(node_coord)

        if i_max != o_max:
            r_dist = (node_coord - i_max) / (o_max-i_max)
            r_dist[r_dist < 0] = 0.0
            rp_dist = 1.0-(node_coord - i_max) / (o_max-i_max)
            rp_dist[r_dist < 0] = 0.0
        else:
            r_dist = rp_dist = numpy.zeros_like(node_coord)

        l_plus_r = l_dist+r_dist
        lp_plus_rp = lp_dist+rp_dist
        return l_plus_r**sigma_exponent, \
               lp_plus_rp**alpha_exponent, \
               l_plus_r**kappa_exponent

    def coefficients_from_boxes(self, discr,
            inner_bbox, outer_bbox=None,
            sigma_magnitude=None, alpha_magnitude=None,
            kappa_magnitude=None, sigma_exponent=None,
            alpha_exponent=None, kappa_exponent=None, dtype=None):

        from math import sqrt, log

        if outer_bbox is None:
            outer_bbox = discr.mesh.bounding_box()

        if sigma_exponent is None:
            sigma_exponent = 2

        if alpha_exponent is None:
            alpha_exponent = 0

        if kappa_exponent is None:
            kappa_exponent = 0

        if sigma_magnitude is None:
            sigma_magnitude = (1+sigma_exponent)*sqrt(self.C[0,0]/self.rho)*log(1/1e-4)*1/(2*(o_min-i_min))

        if alpha_magnitude is None:
            alpha_magnitude = 0.

        if kappa_magnitude is None:
            kappa_magnitude = 1.


        i_min, i_max = inner_bbox
        o_min, o_max = outer_bbox

        from hedge.tools import make_obj_array

        nodes = discr.nodes
        if dtype is not None:
            nodes = nodes.astype(dtype)

        sigma, alpha, kappa = zip(*[self._construct_scalar_coefficients(
            discr, nodes[:,i],
            i_min[i], i_max[i], o_min[i], o_max[i],
            sigma_exponent, alpha_exponent, kappa_exponent)
            for i in range(discr.dimensions)])

        def conv(f):
            return discr.convert_volume(f, kind=discr.compute_kind,
                    dtype=discr.default_scalar_type)

        return self.PMLCoefficients(
                sigma=conv(sigma_magnitude*make_obj_array(sigma)),
                alpha=conv(alpha_magnitude*make_obj_array(alpha)),
                kappa=conv(kappa_magnitude*make_obj_array(kappa)))

    def coefficients_from_width(self, discr, width,
            sigma_magnitude=None, alpha_magnitude=None, kappa_magnitude=None,
            sigma_exponent=None, alpha_exponent=None, kappa_exponent=None,
            dtype=None):
        o_min, o_max = discr.mesh.bounding_box()
        return self.coefficients_from_boxes(discr,
                (o_min+width, o_max-width),
                (o_min, o_max),
                sigma_magnitude, alpha_magnitude, kappa_magnitude,
                sigma_exponent, alpha_exponent, kappa_exponent, dtype)

