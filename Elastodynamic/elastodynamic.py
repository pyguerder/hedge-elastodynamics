# -*- coding: utf-8 -*-
"""Operators for non linear and non-linear Elastodynamics equations."""

from __future__ import division

__authors__ = ["Olivier Bou Matar <olivier.boumatar@iemn.univ-lille1.fr>",
               "Pierre-Yves Guerder <pierre-yves.guerder@centraliens-lille.org>"]
__copyright__ = "Copyright (C) 2010-2011 the authors"
__license__ = "GNU GPLv3 (or more recent equivalent)"


import numpy
from hedge.models import HyperbolicOperator
from pytools import Record
from libraries.utils import Utils
from pymbolic.primitives import IfPositive

def Evaluate(material, val_2, val_1, val_0):
    return IfPositive(material-1,
               val_2,
               IfPositive(material,
                          val_1,
                          val_0))

class ElastoDynamicsOperator(HyperbolicOperator):
    """
    An nD linear Elastodynamics operator.

    dq/dt - dF/dx - dG/dy - dH/dz = 0

    where e.g. in 3D

    q = (rho_v_1, rho_v_2, rho_v_3, F_11, F_22, F_33, 2*F_23, 2*F_13, 2*F_12)
    F = (P_11, P_12, P_13, v_1, 0, 0, 0, v_3, v_2)
    G = (P_12, P_22, P_23, 0, v_2, 0, v_3, 0, v_1)
    H = (P_13, P_23, P_33, 0, 0, v_3, v_2, v_1, 0)

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
                      'open' : 'open' },
                 materials = None,
                 sources = None,
                 flux_type = "lf",
                 ):
        """
        @param sources: should be a table of functions
        that implement
        class:`hedge.data.IFieldDependentGivenFunction`
        or be None.
        @param materials: should be a list
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
        self.sources = sources
        self.flux_type = flux_type
        self.state_null = make_tdep_constant(0)

    def rho_v(self, q):
        return q[0:self.dimensions]

    def F(self, q):
        dim = self.dimensions
        return q[dim:dim+self.dimF[dim]]

    def q(self, w):
        dim = self.dimensions
        return w[1:dim+self.dimF[dim]+1]

    def mat(self, w):
        dim = self.dimensions
        return w[dim+self.dimF[dim]+1]

    def v(self, w):
        from pytools.obj_array import make_obj_array
        q   = self.q(w)
        mat = self.mat(w)
        rho = Evaluate(mat,
                       self.materials[2].rho,
                       self.materials[1].rho,
                       self.materials[0].rho)
        return make_obj_array([rho_v_i/rho for rho_v_i in self.rho_v(q)])

    def P(self, w):
        dim = self.dimensions
        q   = self.q(w)
        mat = self.mat(w)
        Pi = numpy.zeros((self.dimF[dim]), dtype=object)
        F = self.F(q)
        for i in range(self.dimF[dim]):
            for j in range(self.dimF[dim]):
                Ceqij = Evaluate(mat,
                                 self.materials[2].Ceq[i,j],
                                 self.materials[1].Ceq[i,j],
                                 self.materials[0].Ceq[i,j])
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
        fluxes_ph = [fvph[i*n+1:(i+1)*n+1] for i in range(1,d+1)]

        normal = make_normal(d)
        penalty = flux_max(wave_speed_ph.int,wave_speed_ph.ext)*(state_ph.ext-state_ph.int)

        flux_strong = 0.5*sum(n_i*(f_i.ext-f_i.int) for n_i, f_i in zip(normal, fluxes_ph))

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
                        P[0],  # flux rho_v
                        v[0]   # flux F
                        ), "x_flux")
                    ]

        elif dim == 2:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        P[0],P[2],       # flux rho_v
                        v[0],v_null,v[1] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        P[2],P[1],       # flux rho_v
                        v_null,v[1],v[0] # flux F
                        ), "y_flux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        P[0],P[5],P[4],                     # flux rho_v
                        v[0],v_null,v_null,v_null,v[2],v[1] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        P[5],P[1],P[3],                     # flux rho_v
                        v_null,v[1],v_null,v[2],v_null,v[0] # flux F
                        ), "y_flux"),
                    cse(join_fields(
                        P[4],P[3],P[2],                     # flux rho_v
                        v_null,v_null,v[2],v[1],v[0],v_null # flux F
                        ), "z_flux")
                    ]
        else:
            raise ValueError("Invalid dimension")

    def bdry_flux(self, q_bdry, q_null, tag):
        from hedge.tools.symbolic import make_common_subexpression as cse
        from pytools.obj_array import join_fields

        # stress free BCs -------------------------------------------------------
        if tag == self.boundaryconditions_tag['stressfree']:
            signP = -1
            signv = 1
        # fixed BCs -------------------------------------------------------------
        elif tag == self.boundaryconditions_tag['fixed']:
            signP = 1
            signv = -1
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

    def add_sources(self, result):
        from hedge.optemplate import Field
        dim = self.dimensions
        if self.sources is not None:
            if dim == 1:
                result[0] += Field('source_x')
            elif dim == 2:
                result[0] += Field('source_x')
                result[1] += Field('source_y')
            elif dim == 3:
                result[0] += Field('source_x')
                result[1] += Field('source_y')
                result[2] += Field('source_z')
        return result

    def bind_sources(self, t, discr):
        kwargs = {}
        for source in self.sources:
            kwargs[source] = self.sources[source].volume_interpolant(t, discr)
        return kwargs

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

        mat = self.mat(w)
        C00 = Evaluate(mat,
                       self.materials[2].Ceq[0,0],
                       self.materials[1].Ceq[0,0],
                       self.materials[0].Ceq[0,0])
        rho = Evaluate(mat,
                       self.materials[2].rho,
                       self.materials[1].rho,
                       self.materials[0].rho)
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
        nabla = make_nabla(dim)

        result = (numpy.dot(nabla,fluxes) + InverseMassOperator() * (self.flux_num(speed,q,fluxes,bdry_tags_state_and_fluxes)))

        result = self.add_sources(result)

        return result

    def bind(self, discr):
        from hedge.mesh import check_bc_coverage
        check_bc_coverage(discr.mesh, self.boundaryconditions_tag.values())

        compiled_op_template = discr.compile(self.op_template())

        def rhs(t, q):
            extra_kwargs = self.bind_sources(t, discr)
            extra_kwargs['state_null'] = self.state_null.volume_interpolant(t, discr)

            return compiled_op_template(q=q, material=self.material.volume_interpolant(t, discr), **extra_kwargs)

        return rhs

    def max_eigenvalue(self, t, fields=None, discr=None):
        return self.speed


class NLElastoDynamicsOperator(ElastoDynamicsOperator):
    """An nD non linear Elastodynamics operator.

    see YiFeng LI PhD p. 41

    dq/dt - dF/dx - dG/dy - dH/dz = 0

    where e.g. in 3D

    q = (rho_v_1, rho_v_2, rho_v_3, F_11, F_22, F_33, F_23, F_13, F_12, F_32, F_31, F_21)
    F = (P_11, P_21, P_31, v_1, 0, 0, 0, 0, 0, 0, v_3, v_2)
    G = (P_12, P_22, P_32, 0, v_2, 0, 0, 0, v_1, v_3, 0, 0)
    H = (P_13, P_23, P_33, 0, 0, v_3, v_2, v_1, 0, 0, 0, 0)

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
                 sources = None,
                 materials = None,
                 flux_type = "lf",
                 ):
        """
        @param sources: should be a table of functions
        that implement
        class:`hedge.data.IFieldDependentGivenFunction`
        or be None.
        @param materials: should be a list
        of instances of libraries.materials.Material
        """
        from hedge.data import make_tdep_constant

        self.dimensions = dimensions
        self.dimF = [0, 1, 4, 9]

        self.material = material            # The function that gives the material repartition
        self.materials = materials          # The list of used materials
        for i in range(len(materials)):
            self.materials[i].Ceq = Utils.convert_dim(self.materials[i].C, self.dimensions)
        self.speed = speed

        self.boundaryconditions_tag = boundaryconditions_tag
        self.sources = sources
        self.flux_type = flux_type

        self.nonlinearity_type = nonlinearity_type
        self.state_null = make_tdep_constant(0)

    def P(self, w):
        dim = self.dimensions
        q   = w[1:dim+self.dimF[dim]+1]
        Pi = numpy.zeros(self.dimF[dim], dtype=object)
        Ceq = self.C_eq(w)
        F = self.F(q)
        for i in range(self.dimF[dim]):
            for j in range(self.dimF[dim]):
                Pi[i] += Ceq[i,j] * F[j]
        return Pi

    def M(self, w):
        dim = self.dimensions
        mat = w[dim+self.dimF[dim]+1]
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
                                Cnl = Evaluate(mat,
                                               self.materials[2].Cnl[I2,J2,K2],
                                               self.materials[1].Cnl[I2,J2,K2],
                                               self.materials[0].Cnl[I2,J2,K2])
                                Ci2lm = Evaluate(mat,
                                                 self.materials[2].C[I2, LM],
                                                 self.materials[1].C[I2, LM],
                                                 self.materials[0].C[I2, LM])
                                Cilk2 = Evaluate(mat,
                                                 self.materials[2].C[IL, K2],
                                                 self.materials[1].C[IL, K2],
                                                 self.materials[0].C[IL, K2])
                                Ciklm = Evaluate(mat,
                                                 self.materials[2].C[IK, LM],
                                                 self.materials[1].C[IK, LM],
                                                 self.materials[0].C[IK, LM])
                                M[I1,J1,K1] = Cnl \
                                    + Ci2lm * Utils.kronecker(k,n) \
                                    + Cilk2 * Utils.kronecker(j,k) \
                                    + Ciklm * Utils.kronecker(j,n)
        return M

    def C_eq(self, w):
        dim = self.dimensions
        q   = w[1:dim+self.dimF[dim]+1]
        mat = w[dim+self.dimF[dim]+1]
        M = self.M(w)
        if self.nonlinearity_type == "classical":
            C = numpy.zeros((self.dimF[dim],self.dimF[dim]), dtype=object)
            F = self.F(q)
            for i in range(self.dimF[dim]):
                for j in range(self.dimF[dim]):
                    C[i,j] = Evaluate(mat,
                                      self.materials[2].Ceq[i,j],
                                      self.materials[1].Ceq[i,j],
                                      self.materials[0].Ceq[i,j])
                    for k in range(self.dimF[dim]):
                        C[i,j] += 0.5 * M[i, j, k] * F[k]
            return C
        else:
            C = numpy.zeros((self.dimF[dim],self.dimF[dim]), dtype=object)
            F = self.F(q)
            for i in range(self.dimF[dim]):
                for j in range(self.dimF[dim]):
                    Ceqij = Evaluate(mat,
                                     self.materials[2].Ceq[i,j],
                                     self.materials[1].Ceq[i,j],
                                     self.materials[0].Ceq[i,j])
                    for k in range(self.dimF[dim]):
                        C[i,j] = Ceqij
            return C

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
                        P[0],  # flux rho_v
                        v[0]   # flux F
                        ), "x_flux")
                    ]

        elif dim == 2:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        P[0],P[3],              # flux rho_v
                        v[0],v_null,v_null,v[1] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        P[2],P[1],              # flux rho_v
                        v_null,v[1],v[0],v_null # flux F
                        ), "y_flux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        P[0],P[8],P[7],                                          # flux rho_v
                        v[0],v_null,v_null,v_null,v_null,v_null,v_null,v[2],v[1] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        P[5],P[1],P[6],                                          # flux rho_v
                        v_null,v[1],v_null,v_null,v_null,v[0],v[2],v_null,v_null # flux F
                        ), "y_flux"),
                    cse(join_fields(
                        P[4],P[3],P[2],                                          # flux rho_v
                        v_null,v_null,v[2],v[1],v[0],v_null,v_null,v_null,v_null # flux F
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
            signP = -1
            signv = 1
        # fixed BCs -------------------------------------------------------------
        elif tag == self.boundaryconditions_tag['fixed']:
            signP = 1
            signv = -1
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
        self.dimF = [0, 1, 4, 9]

    def F2(self, w):
        dim = self.dimensions
        return w[dim+self.dimF[dim]+2:dim+self.dimF[dim]+2+dim*dim*2]

    def flux(self, w, k):
        from hedge.optemplate import Field
        from hedge.tools.symbolic import make_common_subexpression as cse
        from pytools.obj_array import join_fields

        F2 = self.F2(w) #get F'' from state vector w

        P = self.P(w)
        v = self.v(w)
        v_null = Field('state_null')

        dim = self.dimensions

        if dim == 1:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        (P[0]+F2[0])/k[0],  # flux rho_v
                        (v[0]+F2[1])/k[0]   # flux F
                        ), "x_flux")
                    ]

        elif dim == 2:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        (P[0]+F2[0])/k[0],(P[2]+F2[1])/k[0],       # flux rho_v
                        (v[0]+F2[2])/k[0],v_null,v_null,(v[1]+F2[3])/k[0] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        (P[2]+F2[4])/k[1],(P[1]+F2[5])/k[1],       # flux rho_v
                        v_null,(v[1]+F2[6])/k[1],(v[0]+F2[7])/k[1],v_null # flux F
                        ), "y_flux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        (P[0]+F2[0])/k[0],(P[5]+F2[1])/k[0],(P[4]+F2[2])/k[0],                                             # flux rho_v
                        (v[0]+F2[3])/k[0],v_null,v_null,v_null,v_null,v_null,v_null,(v[2]+F2[4])/k[0],(v[1]+F2[5])/k[0]    # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        (P[5]+F2[6])/k[1],(P[1]+F2[7])/k[1],(P[3]+F2[8])/k[1],                                             # flux rho_v
                        v_null,(v[1]+F2[9])/k[1],v_null,v_null,v_null,(v[2]+F2[10])/k[1],(v[0]+F2[11])/k[1],v_null,v_null  # flux F
                        ), "y_flux"),
                    cse(join_fields(
                        (P[4]+F2[12])/k[2],(P[3]+F2[13])/k[2],(P[2]+F2[14])/k[2],                                          # flux rho_v
                        v_null,v_null,(v[2]+F2[15])/k[2],(v[1]+F2[16])/k[2],(v[0]+F2[17])/k[2],v_null,v_null,v_null,v_null # flux F
                        ), "z_flux")
                    ]
        else:
            raise ValueError("Invalid dimension")


    def bdry_flux(self, q_bdry, q_null, tag):
        from hedge.tools.symbolic import make_common_subexpression as cse
        from pytools.obj_array import join_fields

        # stress free BCs -------------------------------------------------------
        if tag == self.boundaryconditions_tag['stressfree']:
            signP = -1
            signv = 1
        # fixed BCs -------------------------------------------------------------
        elif tag == self.boundaryconditions_tag['fixed']:
            signP = 1
            signv = -1
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
                        signP*P[0],signP*P[2],               # flux rho_v
                        signv*v[0],v_null,v_null,signv*v[1]  # flux F
                        ), "x_bflux"),
                    cse(join_fields(
                        signP*P[2],signP*P[1],               # flux rho_v
                        v_null,signv*v[1],signv*v[0],v_null  # flux F
                        ), "y_bflux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        signP*P[0],signP*P[5],signP*P[4],                                          # flux rho_v
                        signv*v[0],v_null,v_null,v_null,v_null,v_null,signv*v[2],v_null,signv*v[1] # flux F
                        ), "x_bflux"),
                    cse(join_fields(
                        signP*P[5],signP*P[1],signP*P[3],                                          # flux rho_v
                        v_null,signv*v[1],v_null,v_null,signv*v[2],v_null,v_null,signv*v[0],v_null # flux F
                        ), "y_bflux"),
                    cse(join_fields(
                        signP*P[4],signP*P[3],signP*P[2],                                          # flux rho_v
                        v_null,v_null,signv*v[2],signv*v[1],v_null,signv*v[0],v_null,v_null,v_null # flux F
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
        f2 = make_vector_field('f2', dim*dim*2)
        w = join_fields(speed, q, material, f2)

        mat = self.mat(w)
        C00 = Evaluate(mat,
                       self.materials[2].Ceq[0,0],
                       self.materials[1].Ceq[0,0],
                       self.materials[0].Ceq[0,0])
        rho = Evaluate(mat,
                       self.materials[2].rho,
                       self.materials[1].rho,
                       self.materials[0].rho)
        speed = (C00/rho)**0.5

        dim_subset = (True,) * dim + (False,) * (3-dim)

        def pad_vec(v, subset):
            result = numpy.zeros((3,), dtype=object)
            result[numpy.array(subset, dtype=bool)] = v
            return result

        sigma = pad_vec(make_vector_field('sigma', dim),dim_subset)
        alpha = pad_vec(make_vector_field('alpha', dim),dim_subset)
        kappa = pad_vec(make_vector_field('kappa', dim),dim_subset)

        # fluxes ------------------------------------------------------------
        fluxes = self.flux(w, kappa)

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

        bdry_tags_state_and_fluxes = [(tag, bc, self.bdry_flux(bw, bc_null, tag)) 
                                      for tag, bc, bc_null, bw in all_tags_and_bcs]

        # entire operator -----------------------------------------------------
        nabla = make_nabla(dim)

        res_q = (numpy.dot(nabla,fluxes) + InverseMassOperator() * (self.flux_num(speed,q,fluxes,bdry_tags_state_and_fluxes)))

        res_q = self.add_sources(res_q)

        F2 = self.F2(w)
        F = self.F(q)
        
        P = self.P(w)
        v = self.v(w)
        if dim == 1:
            F = [P[0],v[0]]
        elif dim == 2:
            F = [P[0],P[2],v[0],v[1],
                 P[2],P[1],v[1],v[0]]
        elif dim == 3:
            F = [P[0],P[5],P[4],v[0],v[2],v[1],
                 P[5],P[1],P[3],v[1],v[2],v[0],
                 P[4],P[3],P[2],v[2],v[1],v[0]]
        else:
            raise ValueError("Invalid dimension")

        res_f2 = numpy.zeros((dim*dim*2), dtype=object)
        for i in range(dim):
            for j in range(dim*2):
                res_f2[i*dim*2+j] = (-1)*F2[i*dim*2+j]*alpha[i]-sigma[i]/kappa[i]*(F2[i*dim*2+j]+F[i*dim*2+j])

        return join_fields(res_q, res_f2)

    def bind(self, discr, coefs):
        from hedge.mesh import check_bc_coverage
        check_bc_coverage(discr.mesh, self.boundaryconditions_tag.values())

        compiled_op_template = discr.compile(self.op_template())

        def rhs(t, q):
            extra_kwargs = self.bind_sources(t, discr)
            extra_kwargs['state_null'] = self.state_null.volume_interpolant(t, discr)

            dim = self.dimensions            
            return compiled_op_template(q=q[0:dim+self.dimF[dim]],
                                        f2=q[dim+self.dimF[dim]:dim+self.dimF[dim]+dim*dim*2],
                                        material=self.material.volume_interpolant(t, discr),
                                        sigma=coefs.sigma,
                                        alpha=coefs.alpha,
                                        kappa=coefs.kappa,
                                        **extra_kwargs)

        return rhs

    # sigma business ----------------------------------------------------------
    def construct_scalar_coefficients(self, discr, node_coord,
            i_min, i_max, o_min, o_max, sigma_exponent, alpha_exponent, kappa_exponent):
        assert o_min <= i_min <= i_max <= o_max

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

        return l_dist**sigma_exponent, \
               r_dist**sigma_exponent, \
               lp_dist**alpha_exponent, \
               rp_dist**alpha_exponent, \
               l_dist**kappa_exponent, \
               r_dist**kappa_exponent


    def coefficients_from_boxes(self, discr, global_mesh, material,
            inner_bbox, outer_bbox=None,
            sigma_magnitude=None, alpha_magnitude=None,
            kappa_magnitude=None, sigma_exponent=None,
            alpha_exponent=None, kappa_exponent=None, dtype=None):

        from math import sqrt, log

        if outer_bbox is None:
            outer_bbox = global_mesh.bounding_box()

        if sigma_exponent is None:
            sigma_exponent = 2

        if alpha_exponent is None:
            alpha_exponent = 1

        if kappa_exponent is None:
            kappa_exponent = 2

        i_min, i_max = inner_bbox
        o_min, o_max = outer_bbox

        dim = self.dimensions
        sigma_magnitude_l = numpy.zeros(dim)
        sigma_magnitude_r = numpy.zeros(dim)

        for i in range(dim):
            i_l = inner_bbox[0][i]
            i_r = inner_bbox[1][i]
            o_l = outer_bbox[0][i]
            o_r = outer_bbox[1][i]
            sigma_magnitude_l[i] = sigma_magnitude
            if sigma_magnitude is None:
                if i_l-o_l > 0:
                    sigma_magnitude_l[i] = (1+sigma_exponent)*sqrt(material.C[0,0]/material.rho)*log(1/1e-4)*1/(2*(i_l-o_l))
                else:
                    sigma_magnitude_l[i] = 0
    
            sigma_magnitude_r[i] = sigma_magnitude
            if sigma_magnitude is None:
                if o_r-i_r > 0:
                    sigma_magnitude_r[i] = (1+sigma_exponent)*sqrt(material.C[0,0]/material.rho)*log(1/1e-4)*1/(2*(o_r-i_r))
                else:
                    sigma_magnitude_r[i] = 0

        if alpha_magnitude is None:
            alpha_magnitude = 0.

        if kappa_magnitude is None:
            kappa_magnitude = 0.

        from hedge.tools import make_obj_array

        nodes = discr.nodes
        if dtype is not None:
            nodes = nodes.astype(dtype)

        sigma_l_coef, sigma_r_coef, \
        alpha_l_coef, alpha_r_coef, \
        kappa_l_coef, kappa_r_coef = zip(*[self.construct_scalar_coefficients(
            discr, nodes[:,i],
            i_min[i], i_max[i], o_min[i], o_max[i],
            sigma_exponent, alpha_exponent, kappa_exponent)
            for i in range(self.dimensions)])

        def conv(f):
            return discr.convert_volume(f, kind=discr.compute_kind,
                    dtype=discr.default_scalar_type)

        sigma=conv(sigma_magnitude_l*make_obj_array(sigma_l_coef) \
                 + sigma_magnitude_r*make_obj_array(sigma_r_coef))
        alpha=conv(alpha_magnitude*make_obj_array(alpha_l_coef) \
                 + alpha_magnitude*make_obj_array(alpha_r_coef))
        kappa=conv(1+kappa_magnitude*make_obj_array(kappa_l_coef) \
                 + kappa_magnitude*make_obj_array(kappa_r_coef))

        return self.PMLCoefficients(sigma=sigma, alpha=alpha, kappa=kappa)


    def bounding_box(self, mesh):
        try:
            return self._bounding_box
        except AttributeError:
            self._bounding_box = [numpy.min(mesh.points, axis=0), numpy.max(mesh.points, axis=0)]
            return self._bounding_box

    def coefficients_from_width(self, discr, global_mesh, widths, material,
            sigma_magnitude=None, alpha_magnitude=None, kappa_magnitude=None,
            sigma_exponent=None, alpha_exponent=None, kappa_exponent=None,
            dtype=None):
        """
        :param width: [x_l, y_l, z_l, x_r, y_r, z_r]
        """
        o_limits = self.bounding_box(global_mesh)
        i_limits = [numpy.array([o_limits[0][i] + widths[i] for i in range(self.dimensions)]),
                    numpy.array([o_limits[1][i] - widths[i+3] for i in range(self.dimensions)])]

        return self.coefficients_from_boxes(discr,
                global_mesh,
                material,
                i_limits,
                o_limits,
                sigma_magnitude, alpha_magnitude, kappa_magnitude,
                sigma_exponent, alpha_exponent, kappa_exponent, dtype)

class NLNPMLElastoDynamicsOperator(NLElastoDynamicsOperator, NPMLElastoDynamicsOperator):
    """
    Implements NL NPML, based on the 2 previous classes
    """

    def __init__(self, *args, **kwargs):
        #self.add_decay = kwargs.pop("add_decay", True)
        NLElastoDynamicsOperator.__init__(self, *args, **kwargs)

    def flux(self, w, k):
        from hedge.optemplate import Field
        from hedge.tools.symbolic import make_common_subexpression as cse
        from pytools.obj_array import join_fields

        F2 = self.F2(w) #get F'' from state vector w

        P = self.P(w)
        v = self.v(w)
        v_null = Field('state_null')

        dim = self.dimensions

        if dim == 1:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        (P[0]+F2[0])/k[0],  # flux rho_v
                        (v[0]+F2[1])/k[0]   # flux F
                        ), "x_flux")
                    ]

        elif dim == 2:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        (P[0]+F2[0])/k[0],(P[3]+F2[1])/k[0],       # flux rho_v
                        (v[0]+F2[2])/k[0],v_null,v_null,(v[1]+F2[3])/k[0] # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        (P[2]+F2[4])/k[1],(P[1]+F2[5])/k[1],       # flux rho_v
                        v_null,(v[1]+F2[6])/k[1],(v[0]+F2[7])/k[1],v_null # flux F
                        ), "y_flux")
                    ]
        elif dim == 3:
            return [ # one entry for each flux direction
                    cse(join_fields(
                        (P[0]+F2[0])/k[0],(P[8]+F2[1])/k[0],(P[7]+F2[2])/k[0],                                             # flux rho_v
                        (v[0]+F2[3])/k[0],v_null,v_null,v_null,v_null,v_null,v_null,(v[2]+F2[4])/k[0],(v[1]+F2[5])/k[0]    # flux F
                        ), "x_flux"),
                    cse(join_fields(
                        (P[5]+F2[6])/k[1],(P[1]+F2[7])/k[1],(P[6]+F2[8])/k[1],                                             # flux rho_v
                        v_null,(v[1]+F2[9])/k[1],v_null,v_null,v_null,(v[0]+F2[10])/k[1],(v[2]+F2[11])/k[1],v_null,v_null  # flux F
                        ), "y_flux"),
                    cse(join_fields(
                        (P[4]+F2[12])/k[2],(P[3]+F2[13])/k[2],(P[2]+F2[14])/k[2],                                          # flux rho_v
                        v_null,v_null,(v[2]+F2[15])/k[2],(v[1]+F2[16])/k[2],(v[0]+F2[17])/k[2],v_null,v_null,v_null,v_null # flux F
                        ), "z_flux")
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
        f2 = make_vector_field('f2', dim*dim*2)
        w = join_fields(speed, q, material, f2)

        mat = self.mat(w)
        C00 = Evaluate(mat,
                       self.materials[2].Ceq[0,0],
                       self.materials[1].Ceq[0,0],
                       self.materials[0].Ceq[0,0])
        rho = Evaluate(mat,
                       self.materials[2].rho,
                       self.materials[1].rho,
                       self.materials[0].rho)
        speed = (C00/rho)**0.5

        dim_subset = (True,) * dim + (False,) * (3-dim)

        def pad_vec(v, subset):
            result = numpy.zeros((3,), dtype=object)
            result[numpy.array(subset, dtype=bool)] = v
            return result

        sigma = pad_vec(make_vector_field('sigma', dim),dim_subset)
        alpha = pad_vec(make_vector_field('alpha', dim),dim_subset)
        kappa = pad_vec(make_vector_field('kappa', dim),dim_subset)

        # fluxes ------------------------------------------------------------
        fluxes = self.flux(w, kappa)

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

        bdry_tags_state_and_fluxes = [(tag, bc, self.bdry_flux(bw, bc_null, tag)) 
                                      for tag, bc, bc_null, bw in all_tags_and_bcs]

        # entire operator -----------------------------------------------------
        nabla = make_nabla(dim)

        res_q = (numpy.dot(nabla,fluxes) + InverseMassOperator() * (self.flux_num(speed,q,fluxes,bdry_tags_state_and_fluxes)))

        res_q = self.add_sources(res_q)

        F2 = self.F2(w)
        F = self.F(q)
        
        P = self.P(w)
        v = self.v(w)
        if dim == 1:
            F = [P[0],v[0]]
        elif dim == 2:
            F = [P[0],P[3],v[0],v[1],
                 P[2],P[1],v[1],v[0]]
        elif dim == 3:
            F = [P[0],P[8],P[7],v[0],v[2],v[1],
                 P[5],P[1],P[6],v[1],v[2],v[0],
                 P[4],P[3],P[2],v[2],v[1],v[0]]
        else:
            raise ValueError("Invalid dimension")

        res_f2 = numpy.zeros((dim*dim*2), dtype=object)
        for i in range(dim):
            for j in range(dim*2):
                res_f2[i*dim*2+j] = (-1)*F2[i*dim*2+j]*alpha[i]-sigma[i]/kappa[i]*(F2[i*dim*2+j]+F[i*dim*2+j])

        return join_fields(res_q, res_f2)
