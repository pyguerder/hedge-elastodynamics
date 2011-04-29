# -*- coding: utf-8 -*-
"""Operator for non linear Elastodynamics equations."""

from __future__ import division

__copyright__ = """
Copyright (C) 2010-2011:
* Olivier Bou Matar <olivier.boumatar@iemn.univ-lille1.fr>
* Pierre-Yves Guerder <pierre-yves.guerder@centraliens-lille.org>
"""

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
from hedge.models import HyperbolicOperator
from pytools import Record
from libraries.utils import Utils


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
                 material,
                 speed,
                 nonlinearity_type = "classical",
                 boundaryconditions_tag = \
                    { 'stressfree' : 'stressfree',
                      'fixed' : 'fixed',
                      'open' : 'open' },
                 source = None,
                 materials = None,
                 flux_type = "lf",
                 ):
        """
        :param source: should implement
        class:`hedge.data.IFieldDependentGivenFunction`
        or be None.
        :param materials: should be a list
        of instances of libraries.materials.Material
        """
        from hedge.data import make_tdep_constant

        self.dimensions = dimensions

        self.material = material            # The function that gives the material repartition
        self.materials = materials          # The list of used materials
        for i in range(len(materials)):
            self.materials[i].Ceq = Utils.convert_dim(self.materials[i].C, self.dimensions)
        self.speed = speed

        self.boundaryconditions_tag = boundaryconditions_tag
        self.source = source
        self.flux_type = flux_type

        self.nonlinearity_type = nonlinearity_type
        self.state_null = make_tdep_constant(0)

    def rho_v(self, q):
        return q[0:self.dimensions]

    def F(self, q):
        return q[self.dimensions:(self.dimensions+1)*self.dimensions]

    def v(self, w):
        from pytools.obj_array import make_obj_array
        from pymbolic.primitives import IfPositive
        q   = w[1:self.dimensions*(self.dimensions+1)+1]
        mat = w[self.dimensions*(self.dimensions+1)+1]
        rho = IfPositive(mat, self.materials[0].rho, self.materials[1].rho)
        return make_obj_array([rho_v_i/rho for rho_v_i in self.rho_v(q)])

    def P(self, w):
        q   = w[1:self.dimensions*(self.dimensions+1)+1]
        Pi = numpy.zeros(self.dimensions**2, dtype=object)
        Ceq = self.C_eq(w)
        F = self.F(q)
        for i in range(self.dimensions**2):
            for j in range(self.dimensions**2):
                Pi[i] += Ceq[i,j] * F[j]
        return Pi

    def M(self, w):
        dim = self.dimensions
        from pymbolic.primitives import IfPositive
        mat = w[self.dimensions*(self.dimensions+1)+1]
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
                                Cnl = IfPositive(mat, self.materials[0].Cnl[I2,J2,K2], self.materials[1].Cnl[I2,J2,K2])
                                Ci2lm = IfPositive(mat, self.materials[0].C[I2, LM], self.materials[1].C[I2, LM])
                                Cilk2 = IfPositive(mat, self.materials[0].C[IL, K2], self.materials[1].C[IL, K2])
                                Ciklm = IfPositive(mat, self.materials[0].C[IK, LM], self.materials[1].C[IK, LM])
                                M[I1,J1,K1] = Cnl \
                                    + Ci2lm * Utils.kronecker(k,n) \
                                    + Cilk2 * Utils.kronecker(j,k) \
                                    + Ciklm * Utils.kronecker(j,n)
        return M

    def C_eq(self, w):
        from pymbolic.primitives import IfPositive
        q   = w[1:self.dimensions*(self.dimensions+1)+1]
        mat = w[self.dimensions*(self.dimensions+1)+1]
        M = self.M(w)
        if self.nonlinearity_type == "classical":
            C = numpy.zeros((self.dimensions**2,self.dimensions**2), dtype=object)
            F = self.F(q)
            for i in range(self.dimensions**2):
                for j in range(self.dimensions**2):
                    Ceqij = IfPositive(mat, self.materials[0].Ceq[i,j], self.materials[1].Ceq[i,j])
                    for k in range(self.dimensions**2):
                        C[i,j] += Ceqij + 0.5 * M[i, j, k] * F[k]
            return C
        else:
            Ceq = IfPositive(mat, self.materials[0].Ceq, self.materials[1].Ceq)
            return Ceq

    def flux_num(self, wave_speed, q, fluxes, bdry_tags_state_and_fluxes):
        from hedge.flux import FluxVectorPlaceholder, make_normal, flux_max
        from hedge.optemplate import get_flux_operator, BoundaryPair
        from pytools.obj_array import join_fields

        n = len(q)
        d = len(fluxes)
        fvph = FluxVectorPlaceholder(n*(d+1)+1)
        wave_speed_ph = fvph[0]
        state_ph  = fvph[1:n+1]
        # ??? modifié le range(1,d+1) en range(1,d) ; sinon : index out of bound
        fluxes_ph = [fvph[i*n+1:(i+1)*n+1] for i in range(1,d)]

        normal = make_normal(d)
        penalty = flux_max(wave_speed_ph.int,wave_speed_ph.ext)*(state_ph.ext-state_ph.int)

        flux_strong = 0.5*sum(n_i*(f_i.int-f_i.ext) for n_i, f_i in zip(normal, fluxes_ph))

        if self.flux_type == "central":
            pass
        elif self.flux_type == "lf":
            flux_strong = 0.5*penalty + flux_strong
        else:
            raise ValueError("Invalid flux type '%s'" % self.flux_type)

        flux_op = get_flux_operator(flux_strong)
        int_operand = join_fields(wave_speed,q,*fluxes)

        return (flux_op(int_operand)
                +sum(flux_op(BoundaryPair(int_operand, join_fields(0, bdry_state, *bdry_fluxes), tag))
                     for tag, bdry_state, bdry_fluxes in bdry_tags_state_and_fluxes))

    def flux(self, q):
        from hedge.optemplate import Field
        from hedge.tools.symbolic import make_common_subexpression as cse
        from pytools.obj_array import join_fields

        P = self.P(q)
        v = self.v(q)
        v_null = Field('state_null')

        dim = self.dimensions

        if dim == 1:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        -P[0],  # flux rho_v
                        -v[0]   # flux F
                        ), "x_flux")
                    ]

        elif dim == 2:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        -P[0],-P[3],              # flux rho_v
                        -v[0],v_null,v_null,-v[1] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        -P[2],-P[1],              # flux rho_v
                        v_null,-v[1],-v[0],v_null # flux F
                        ), "y_flux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        -P[0],-P[8],-P[7],                                          # flux rho_v
                        -v[0],v_null,v_null,v_null,v_null,v_null,v_null,-v[2],-v[1] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        -P[5],-P[1],-P[6],                                          # flux rho_v
                        v_null,-v[1],v_null,v_null,v_null,-v[0],-v[2],v_null,v_null # flux F
                        ), "y_flux"),
                    cse(join_fields(
                        -P[4],-P[3],-P[2],                                          # flux rho_v
                        v_null,v_null,-v[2],-v[1],-v[0],v_null,v_null,v_null,v_null # flux F
                        ), "z_flux")
                    ]
        else:
            raise ValueError("Invalid dimension")

    def bdry_flux(self, q_bdry, q_null, tag):
        from hedge.tools.symbolic import make_common_subexpression as cse
        from pytools.obj_array import join_fields

        dim = self.dimensions

        # stress free BCs -------------------------------------------------------
        if tag == self.boundaryconditions_tag['stressfree']:
            signP = 1
            signv = -1
        # fixed BCs -------------------------------------------------------------
        elif tag == self.boundaryconditions_tag['fixed']:
            signP = -1
            signv = 1
        else:
            raise ValueError("Invalid boundary conditions")

        P = self.P(q_bdry)
        v = self.v(q_bdry)
        v_null = q_null

        if dim == 1:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        signP*P[0], # flux rho_v
                        signv*v[0]  # flux F
                        ), "x_bflux")
                    ]
        elif dim == 2:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        signP*P[0],signP*P[3],              # flux rho_v
                        signv*v[0],v_null,v_null,signv*v[1] # flux F
                        ), "x_bflux"),
                    cse(join_fields(
                        signP*P[2],signP*P[1],              # flux rho_v
                        v_null,signv*v[1],signv*v[0],v_null # flux F
                        ), "y_bflux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        signP*P[0],signP*P[8],signP*P[7],                                          # flux rho_v
                        signv*v[0],v_null,v_null,v_null,v_null,v_null,v_null,signv*v[2],signv*v[1] # flux F
                        ), "x_bflux"),
                    cse(join_fields(
                        signP*P[5],signP*P[1],signP*P[6],                                          # flux rho_v
                        v_null,signv*v[1],v_null,v_null,v_null,signv*v[0],signv*v[2],v_null,v_null # flux F
                        ), "y_bflux"),
                    cse(join_fields(
                        signP*P[4],signP*P[3],signP*P[2],                                          # flux rho_v
                        v_null,v_null,signv*v[2],signv*v[1],signv*v[0],v_null,v_null,v_null,v_null # flux F
                        ), "z_bflux")
                    ]
        else:
            raise ValueError("Invalid dimension")

    def op_template(self):
        from hedge.optemplate import \
                Field, \
                make_nabla, \
                InverseMassOperator, \
                BoundarizeOperator
        from hedge.optemplate.tools import make_vector_field
        from pytools.obj_array import join_fields

        dim = self.dimensions

        speed = Field('speed')
        q = make_vector_field('q', dim*(dim+1))
        material = Field('material')
        w = join_fields(speed,q,material)

        from pymbolic.primitives import IfPositive
        mat = w[self.dimensions*(self.dimensions+1)]
        C00 = IfPositive(mat, self.materials[0].Ceq[0,0], self.materials[1].Ceq[0,0])
        rho = IfPositive(mat, self.materials[0].rho, self.materials[1].rho)
        speed = (C00/rho)**0.5


        # fluxes ------------------------------------------------------------
        fluxes = self.flux(w)


        # boundary conditions ---------------------------------------------------
        state_bc_stressfree = BoundarizeOperator(self.boundaryconditions_tag['stressfree'])(q)
        state_bc_stressfree_null = BoundarizeOperator(self.boundaryconditions_tag['stressfree'])(0)
        state_bc_fixed = BoundarizeOperator(self.boundaryconditions_tag['fixed'])(q)
        state_bc_fixed_null = BoundarizeOperator(self.boundaryconditions_tag['fixed'])(0)
        w_bc_stressfree = BoundarizeOperator(self.boundaryconditions_tag['stressfree'])(w)
        w_bc_fixed = BoundarizeOperator(self.boundaryconditions_tag['fixed'])(w)

        all_tags_and_bcs = [
                (self.boundaryconditions_tag['stressfree'], state_bc_stressfree, state_bc_stressfree_null, w_bc_stressfree),
                (self.boundaryconditions_tag['fixed'], state_bc_fixed, state_bc_fixed_null, w_bc_fixed)
                           ]

        bdry_tags_state_and_fluxes = [(tag, bc, self.bdry_flux(bw, bc_null, tag)) for tag, bc, bc_null, bw in all_tags_and_bcs]

        # entire operator -----------------------------------------------------
        from math import sin, cos, pi
        nabla = make_nabla(dim)

        result = (-numpy.dot(nabla,fluxes) + InverseMassOperator() * (self.flux_num(speed,q,fluxes,bdry_tags_state_and_fluxes)))

        if self.source is not None:
            result[0] += Field('source_v_x')*sin(10*pi/180)
            result[1] += Field('source_v_x')*cos(10*pi/180)

        return result

    def bind(self, discr):
        from hedge.mesh import check_bc_coverage
        check_bc_coverage(discr.mesh, self.boundaryconditions_tag.values())

        compiled_op_template = discr.compile(self.op_template())

        def rhs(t, q):
            extra_kwargs = {}

            if self.source is not None:
                extra_kwargs['source_v_x'] = self.source.volume_interpolant(t, discr)
            extra_kwargs['state_null'] = self.state_null.volume_interpolant(t, discr)

            return compiled_op_template(q=q, material=self.material.volume_interpolant(t, discr), **extra_kwargs)

        return rhs

    def max_eigenvalue(self, t, fields=None, discr=None):
        return self.speed


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
    def __init__(self,
                 dimensions,
                 material,
                 speed,
                 boundaryconditions_tag = \
                    { 'stressfree' : 'stressfree',
                      'fixed' : 'fixed',
                      'open' : 'open',
                      'materials' : None },
                 materials = None,
                 source = None,
                 flux_type = "lf",
                 ):
        """
        :param source: should implement
        class:`hedge.data.IFieldDependentGivenFunction`
        or be None.
        :param materials: should be a list
        of instances of libraries.materials.Material
        """
        from hedge.data import make_tdep_constant

        self.dimensions = dimensions
        self.dimF = [0, 1, 3, 6]

        self.material = material            # The function that gives the material repartition
        self.materials = materials          # The list of used materials
        for i in range(len(materials)):
            self.materials[i].Ceq = Utils.convert_dim(self.materials[i].C, self.dimensions)
        self.speed = speed

        self.boundaryconditions_tag = boundaryconditions_tag
        self.source = source
        self.flux_type = flux_type
        self.state_null = make_tdep_constant(0)

    def rho_v(self, q):
        return q[0:self.dimensions]

    def F(self, q):
        return q[self.dimensions:self.dimensions+self.dimF[self.dimensions]]

    def v(self, w):
        from pytools.obj_array import make_obj_array
        from pymbolic.primitives import IfPositive
        q   = w[1:self.dimensions+self.dimF[self.dimensions]+1]
        mat = w[self.dimensions+self.dimF[self.dimensions]+1]
        rho = IfPositive(mat, self.materials[0].rho, self.materials[1].rho)
        return make_obj_array([rho_v_i/rho for rho_v_i in self.rho_v(q)])

    def P(self, w):
        from pymbolic.primitives import IfPositive
        q   = w[1:self.dimensions+self.dimF[self.dimensions]+1]
        mat = w[self.dimensions+self.dimF[self.dimensions]+1]
        Pi = numpy.zeros((self.dimF[self.dimensions]), dtype=object)
        F = self.F(q)
        for i in range(self.dimF[self.dimensions]):
            for j in range(self.dimF[self.dimensions]):
                Ceqij = IfPositive(mat, self.materials[0].Ceq[i,j], self.materials[1].Ceq[i,j])
                Pi[i] += Ceqij*F[j]
        return Pi

    def flux_num(self, wave_speed, q, fluxes, bdry_tags_state_and_fluxes):
        from hedge.flux import FluxVectorPlaceholder, make_normal, flux_max
        from hedge.optemplate import get_flux_operator, BoundaryPair
        from pytools.obj_array import join_fields

        n = len(q)
        d = len(fluxes)
        fvph = FluxVectorPlaceholder(n*(d+1)+1)
        wave_speed_ph = fvph[0]
        state_ph  = fvph[1:n+1]
        # ??? range(1,d+1) modifié en range(1,d) ; sinon : index out of bound
        fluxes_ph = [fvph[i*n+1:(i+1)*n+1] for i in range(1,d+1)]

        normal = make_normal(d)
        penalty = flux_max(wave_speed_ph.int,wave_speed_ph.ext)*(state_ph.ext-state_ph.int)

        flux_strong = 0.5*sum(n_i*(f_i.int-f_i.ext) for n_i, f_i in zip(normal, fluxes_ph))

        if self.flux_type == "central":
            pass
        elif self.flux_type == "lf":
            flux_strong = 0.5*penalty + flux_strong
        else:
            raise ValueError("Invalid flux type '%s'" % self.flux_type)

        flux_op = get_flux_operator(flux_strong)
        int_operand = join_fields(wave_speed,q,*fluxes)

        return (flux_op(int_operand)
                +sum(flux_op(BoundaryPair(int_operand, join_fields(0,bdry_state, *bdry_fluxes), tag))
                     for tag, bdry_state, bdry_fluxes in bdry_tags_state_and_fluxes))


    def flux(self, w):
        from hedge.optemplate import Field
        from hedge.tools.symbolic import make_common_subexpression as cse
        from pytools.obj_array import join_fields

        P = self.P(w)
        v = self.v(w)
        v_null = Field('state_null')

        dim = self.dimensions

        if dim == 1:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        -P[0],  # flux rho_v
                        -v[0]   # flux F
                        ), "x_flux")
                    ]

        elif dim == 2:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        -P[0],-P[2],       # flux rho_v
                        -v[0],v_null,-v[1] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        -P[2],-P[1],       # flux rho_v
                        v_null,-v[1],-v[0] # flux F
                        ), "y_flux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        -P[0],-P[5],-P[4],                     # flux rho_v
                        -v[0],v_null,v_null,v_null,-v[2],-v[1] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        -P[5],-P[1],-P[3],                     # flux rho_v
                        v_null,-v[1],v_null,-v[2],v_null,-v[0] # flux F
                        ), "y_flux"),
                    cse(join_fields(
                        -P[4],-P[3],-P[2],                     # flux rho_v
                        v_null,v_null,-v[2],-v[1],-v[0],v_null # flux F
                        ), "z_flux")
                    ]
        else:
            raise ValueError("Invalid dimension")

    def bdry_flux(self, q_bdry, q_null, tag):
        from hedge.tools.symbolic import make_common_subexpression as cse
        from pytools.obj_array import join_fields

        # stress free BCs -------------------------------------------------------
        if tag == self.boundaryconditions_tag['stressfree']:
            signP = 1
            signv = -1
        # fixed BCs -------------------------------------------------------------
        elif tag == self.boundaryconditions_tag['fixed']:
            signP = -1
            signv = 1
        else:
            raise ValueError("Invalid boundary conditions")

        dim = self.dimensions

        P = self.P(q_bdry)
        v = self.v(q_bdry)
        v_null = q_null

        if dim == 1:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        signP*P[0], # flux rho_v
                        signv*v[0]  # flux F
                        ), "x_bflux")
                    ]
        elif dim == 2:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        signP*P[0],signP*P[2],        # flux rho_v
                        signv*v[0],v_null,signv*v[1]  # flux F
                        ), "x_bflux"),
                    cse(join_fields(
                        signP*P[2],signP*P[1],        # flux rho_v
                        v_null,signv*v[1],signv*v[0]  # flux F
                        ), "y_bflux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        signP*P[0],signP*P[5],signP*P[4],                     # flux rho_v
                        signv*v[0],v_null,v_null,v_null,signv*v[2],signv*v[1] # flux F
                        ), "x_bflux"),
                    cse(join_fields(
                        signP*P[5],signP*P[1],signP*P[3],                     # flux rho_v
                        v_null,signv*v[1],v_null,signv*v[2],v_null,signv*v[0] # flux F
                        ), "y_bflux"),
                    cse(join_fields(
                        signP*P[4],signP*P[3],signP*P[2],                     # flux rho_v
                        v_null,v_null,signv*v[2],signv*v[1],signv*v[0],v_null # flux F
                        ), "z_bflux")
                    ]
        else:
            raise ValueError("Invalid dimension")


    def op_template(self):
        from hedge.optemplate import \
                Field, \
                make_nabla, \
                InverseMassOperator, \
                BoundarizeOperator
        from hedge.optemplate.tools import make_vector_field
        from pytools.obj_array import join_fields

        dim = self.dimensions

        speed = self.speed
        q = make_vector_field('q', dim+self.dimF[dim])
        material = Field('material')
        w = join_fields(speed,q,material)

        from pymbolic.primitives import IfPositive
        mat = w[self.dimensions+self.dimF[self.dimensions]+1]
        C00 = IfPositive(mat, self.materials[0].Ceq[0,0], self.materials[1].Ceq[0,0])
        rho = IfPositive(mat, self.materials[0].rho, self.materials[1].rho)
        speed = (C00/rho)**0.5


        # fluxes ------------------------------------------------------------
        fluxes = self.flux(w)


        # boundary conditions ---------------------------------------------------
        q_bc_stressfree = BoundarizeOperator(self.boundaryconditions_tag['stressfree'])(q)
        q_bc_stressfree_null = BoundarizeOperator(self.boundaryconditions_tag['stressfree'])(0)
        q_bc_fixed = BoundarizeOperator(self.boundaryconditions_tag['fixed'])(q)
        q_bc_fixed_null = BoundarizeOperator(self.boundaryconditions_tag['fixed'])(0)
        w_bc_stressfree = BoundarizeOperator(self.boundaryconditions_tag['stressfree'])(w)
        w_bc_fixed = BoundarizeOperator(self.boundaryconditions_tag['fixed'])(w)

        all_tags_and_bcs = [
                (self.boundaryconditions_tag['stressfree'], q_bc_stressfree, q_bc_stressfree_null, w_bc_stressfree),
                (self.boundaryconditions_tag['fixed'], q_bc_fixed, q_bc_fixed_null, w_bc_fixed)
                           ]

        bdry_tags_state_and_fluxes = [(tag, bc, self.bdry_flux(bw, bc_null, tag)) for tag, bc, bc_null, bw in all_tags_and_bcs]

        # entire operator -----------------------------------------------------
        from math import sin, cos, pi
        nabla = make_nabla(dim)

        result = (-numpy.dot(nabla,fluxes) + InverseMassOperator() * (self.flux_num(speed,q,fluxes,bdry_tags_state_and_fluxes)))

        if self.source is not None:
            result[0] += Field('source_v_x')*sin(10*pi/180)
            result[1] += Field('source_v_x')*cos(10*pi/180)

        return result

    def bind(self, discr):
        from hedge.mesh import check_bc_coverage
        check_bc_coverage(discr.mesh, self.boundaryconditions_tag.values())

        compiled_op_template = discr.compile(self.op_template())

        def rhs(t, q):
            extra_kwargs = {}

            if self.source is not None:
                extra_kwargs['source_v_x'] = self.source.volume_interpolant(t, discr)
            extra_kwargs['state_null'] = self.state_null.volume_interpolant(t, discr)

            return compiled_op_template(q=q, material=self.material.volume_interpolant(t, discr), **extra_kwargs)

        return rhs

    def max_eigenvalue(self, t, fields=None, discr=None):
        return self.speed


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

        #from hedge.optemplate import make_vector_field
        #sigma = make_vector_field("sigma", dim)
        #alpha = make_vector_field("alpha", dim)
        #kappa = make_vector_field("kappa", dim)


        rhs = numpy.zeros((2*dim**2,), dtype=object)

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
            q = make_vector_field("q", self.dimensions+self.dimF[self.dimensions])

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

        i_min, i_max = inner_bbox
        o_min, o_max = outer_bbox

        if sigma_magnitude is None:
            sigma_magnitude = (1+sigma_exponent)*sqrt(self.C[0,0]/self.rho)*log(1/1e-4)*1/(2*(o_min-i_min))

        if alpha_magnitude is None:
            alpha_magnitude = 0.

        if kappa_magnitude is None:
            kappa_magnitude = 1.

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

