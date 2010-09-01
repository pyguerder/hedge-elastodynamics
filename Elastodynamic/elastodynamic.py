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
    def __init__(self, dimensions, rho,
            C, stressfree_tag="stressfree",
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
        self.rho = rho

        self.stressfree_tag = stressfree_tag
        self.fixed_tag = fixed_tag
        self.open_tag = open_tag

        self.source = source
	self.flux_type = flux_type

    def rho_v(self, state):
        return state[0:self.dimensions]

    def F(self, state):
        return state[self.dimensions:(self.dimensions+1)*self.dimensions]

    def v(self, state):
        from hedge.tools import make_obj_array
        return make_obj_array([
                rho_v_i/self.rho
                for rho_v_i in self.rho_v(state)])

    def P(self, state):
        Pi = numpy.zeros((self.dimensions**2), dtype=object)
	for i in range(self.dimensions**2):
	    for j in range(self.dimensions**2):
	        Pi[i] += self.C[i,j]*self.F(state)[j]
        return Pi

    def flux_num(self, state, fluxes, bdry_tags_state_and_fluxes):
        from hedge.flux import FluxVectorPlaceholder, make_normal
        from hedge.optemplate import get_flux_operator, BoundaryPair
        from hedge.tools import join_fields
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
                make_vector_field, \
                BoundaryPair, \
                get_flux_operator, \
                make_nabla, \
                InverseMassOperator, \
                BoundarizeOperator

	dim = self.dimensions

        q = make_vector_field("q", dim*(dim+1))
	rho_v = q[0:dim]
	F = q[dim:dim*(dim+1)]

        def flux(q):
	    from hedge.tools.symbolic import make_common_subexpression as cse
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

	def bdry_flux(q_bdry, q_vol, tag):
	    from hedge.tools.symbolic import make_common_subexpression as cse

            # stress free BCs -------------------------------------------------------
	    if tag == self.stressfree_tag:
 	        if dim == 1:
                    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[0],

                                # flux F
                                -self.v(q_bdry)[0]
                                ), "x_flux")
                            ]

	        elif dim == 2:
                    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[0],self.P(q_bdry)[3],

                                # flux F
                                -self.v(q_bdry)[0],0,0,-self.v(q_bdry)[1]
                                ), "x_flux"),
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[2],self.P(q_bdry)[1],

                                # flux F
			        0,-self.v(q_bdry)[1],-self.v(q_bdry)[0],0
                                ), "y_flux")
                            ]
	        else:
		    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[0],self.P(q_bdry)[8],self.P(q_bdry)[7],

                                # flux F
                                -self.v(q_bdry)[0],0,0,0,0,0,0,-self.v(q_bdry)[2],-self.v(q_bdry)[1]
                                ), "x_flux"),
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[5],self.P(q_bdry)[1],self.P(q_bdry)[6],

                                # flux F
			        0,-self.v(q_bdry)[1],0,0,0,-self.v(q_bdry)[0],-self.v(q_bdry)[2],0,0
                                ), "y_flux"),
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[4],self.P(q_bdry)[3],self.P(q_bdry)[2],

                                # flux F
			        0,0,-self.v(q_bdry)[2],-self.v(q_bdry)[1],-self.v(q_bdry)[0],0,0,0,0
                                ), "z_flux"),
                            ]
            # fixed BCs ---------------------------------------------------------
	    elif tag == self.fixed_tag:
 	        if dim == 1:
                    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[0],

                                # flux F
                                self.v(q_bdry)[0]
                                ), "x_flux")
                            ]

	        elif dim == 2:
                    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[0],-self.P(q_bdry)[3],

                                # flux F
                                self.v(q_bdry)[0],0,0,self.v(q_bdry)[1]
                                ), "x_flux"),
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[2],-self.P(q_bdry)[1],

                                # flux F
			        0,self.v(q_bdry)[1],self.v(q_bdry)[0],0
                                ), "y_flux")
                            ]
	        else:
		    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[0],-self.P(q_bdry)[8],-self.P(q_bdry)[7],

                                # flux F
                                self.v(q_bdry)[0],0,0,0,0,0,0,self.v(q_bdry)[2],self.v(q_bdry)[1]
                                ), "x_flux"),
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[5],-self.P(q_bdry)[1],-self.P(q_bdry)[6],

                                # flux F
			        0,self.v(q_bdry)[1],0,0,0,self.v(q_bdry)[0],self.v(q_bdry)[2],0,0
                                ), "y_flux"),
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[4],-self.P(q_bdry)[3],-self.P(q_bdry)[2],

                                # flux F
			        0,0,self.v(q_bdry)[2],self.v(q_bdry)[1],self.v(q_bdry)[0],0,0,0,0
                                ), "z_flux")
                            ]
	    else:
	        raise ValueError("invalid Bdry conditions")

        # fluxes ------------------------------------------------------------
	from hedge.tools import join_fields
 	
	fluxes = flux(q)

        # boundary conditions ---------------------------------------------------
	state_bc_stressfree = BoundarizeOperator(self.stressfree_tag) * q
	state_bc_fixed      = BoundarizeOperator(self.fixed_tag) * q

	all_tags_and_bcs = [
		(self.stressfree_tag, state_bc_stressfree),
		(self.fixed_tag, state_bc_fixed)
			   ]

	bdry_tags_state_and_fluxes = [
		(tag, bc, bdry_flux(bc, q, tag))
		for tag, bc in all_tags_and_bcs]

        # entire operator -----------------------------------------------------
        nabla = make_nabla(dim)

        result = (-numpy.dot(nabla,fluxes)
                  +
                  InverseMassOperator() * (
			self.flux_num(q,fluxes,bdry_tags_state_and_fluxes)
					   ))

        if self.source is not None:
            result[0] += Field("source_v_x")
	    #result[1] += Field("source_v_x")

        return result

    def bind(self, discr):
        from hedge.mesh import check_bc_coverage
        check_bc_coverage(discr.mesh, [
            self.stressfree_tag,
            self.fixed_tag,
            self.open_tag])

        compiled_op_template = discr.compile(self.op_template())

        def rhs(t, q):
            extra_kwargs = {"q": q}

            if self.source is not None:
                extra_kwargs["source_v_x"] = self.source.volume_interpolant(t, discr)

            return compiled_op_template(**extra_kwargs)

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
        self.rho = rho

        self.stressfree_tag = stressfree_tag
        self.fixed_tag = fixed_tag
        self.open_tag = open_tag

        self.source = source
	self.flux_type = flux_type

    def rho_v(self, state):
        return state[0:self.dimensions]

    def F(self, state):
        return state[self.dimensions:self.dimensions+self.dimF[self.dimensions]]

    def v(self, state):
        from hedge.tools import make_obj_array
        return make_obj_array([
                rho_v_i/self.rho
                for rho_v_i in self.rho_v(state)])

    def P(self, state):
        Pi = numpy.zeros((self.dimF[self.dimensions]), dtype=object)
	for i in range(self.dimF[self.dimensions]):
	    for j in range(self.dimF[self.dimensions]):
	            Pi[i] += self.C[i,j]*self.F(state)[j]
        return Pi

    def flux_num(self, state, fluxes, bdry_tags_state_and_fluxes):
        from hedge.flux import FluxVectorPlaceholder, make_normal
        from hedge.optemplate import get_flux_operator, BoundaryPair
        from hedge.tools import join_fields
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
                make_vector_field, \
                BoundaryPair, \
                get_flux_operator, \
                make_nabla, \
                InverseMassOperator, \
                BoundarizeOperator

	dim = self.dimensions

        q = make_vector_field("q", dim+self.dimF[dim])
	rho_v = q[0:dim]
	F = q[dim:dim+self.dimF[dim]]

        def flux(q):
	    from hedge.tools.symbolic import make_common_subexpression as cse
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
                            -self.v(q)[0],0.0,-self.v(q)[1]
                            ), "x_flux"),
                        cse(join_fields(
                            # flux rho_v
                            -self.P(q)[2],-self.P(q)[1],

                            # flux F
			    0.0,-self.v(q)[1],-self.v(q)[0]
                            ), "y_flux")
                        ]
	    else:
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

	def bdry_flux(q_bdry, q_vol, tag):
	    from hedge.tools.symbolic import make_common_subexpression as cse

            # stress free BCs -------------------------------------------------------
	    if tag == self.stressfree_tag:
 	        if dim == 1:
                    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[0],

                                # flux F
                                -self.v(q_bdry)[0]
                                ), "x_flux")
                            ]

	        elif dim == 2:
                    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[0],self.P(q_bdry)[2],

                                # flux F
                                -self.v(q_bdry)[0],0.0,-self.v(q_bdry)[1]
                                ), "x_flux"),
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[2],self.P(q_bdry)[1],

                                # flux F
			        0.0,-self.v(q_bdry)[1],-self.v(q_bdry)[0]
                                ), "y_flux")
                            ]
	        else:
		    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[0],self.P(q_bdry)[5],self.P(q_bdry)[4],

                                # flux F
                                -self.v(q_bdry)[0],0,0,0,-self.v(q_bdry)[2],-self.v(q_bdry)[1]
                                ), "x_flux"),
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[5],self.P(q_bdry)[1],self.P(q_bdry)[3],

                                # flux F
			        0,-self.v(q_bdry)[1],0,-self.v(q_bdry)[2],0,-self.v(q_bdry)[0]
                                ), "y_flux"),
                            cse(join_fields(
                                # flux rho_v
                                self.P(q_bdry)[4],self.P(q_bdry)[3],self.P(q_bdry)[2],

                                # flux F
			        0,0,-self.v(q_bdry)[2],-self.v(q_bdry)[1],-self.v(q_bdry)[0],0
                                ), "z_flux")
                            ]
            # fixed BCs ---------------------------------------------------------
	    elif tag == self.fixed_tag:
 	        if dim == 1:
                    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[0],

                                # flux F
                                self.v(q_bdry)[0]
                                ), "x_flux")
                            for i in range(self.dimensions)]

	        elif dim == 2:
                    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[0],-self.P(q_bdry)[2],

                                # flux F
                                self.v(q_bdry)[0],0.0,self.v(q_bdry)[1]
                                ), "x_flux"),
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[2],-self.P(q_bdry)[1],

                                # flux F
			        0.0,self.v(q_bdry)[1],self.v(q_bdry)[0]
                                ), "y_flux")
                            ]
	        else:
		    return [ # one entry for each flux direction
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[0],-self.P(q_bdry)[5],-self.P(q_bdry)[4],

                                # flux F
                                self.v(q_bdry)[0],0,0,0,self.v(q_bdry)[2],self.v(q_bdry)[1]
                                ), "x_flux"),
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[5],-self.P(q_bdry)[1],-self.P(q_bdry)[3],

                                # flux F
			        0,self.v(q_bdry)[1],0,self.v(q_bdry)[2],0,self.v(q_bdry)[0]
                                ), "y_flux"),
                            cse(join_fields(
                                # flux rho_v
                                -self.P(q_bdry)[4],-self.P(q_bdry)[3],-self.P(q_bdry)[2],

                                # flux F
			        0,0,self.v(q_bdry)[2],self.v(q_bdry)[1],self.v(q_bdry)[0],0
                                ), "z_flux")
                            ]
	    else:
	        raise ValueError("invalid Bdry conditions")

        # fluxes ------------------------------------------------------------
	from hedge.tools import join_fields
 	
	fluxes = flux(q)

        # boundary conditions ---------------------------------------------------
	state_bc_stressfree = BoundarizeOperator(self.stressfree_tag) * q
	state_bc_fixed      = BoundarizeOperator(self.fixed_tag) * q

	all_tags_and_bcs = [
		(self.stressfree_tag, state_bc_stressfree),
		(self.fixed_tag, state_bc_fixed)
			   ]

	bdry_tags_state_and_fluxes = [
		(tag, bc, bdry_flux(bc, q, tag))
		for tag, bc in all_tags_and_bcs]

        # entire operator -----------------------------------------------------
        nabla = make_nabla(dim)

        result = (-numpy.dot(nabla,fluxes)
                  +
                  InverseMassOperator() * (
			self.flux_num(q,fluxes,bdry_tags_state_and_fluxes)
					   ))

        if self.source is not None:
            result[0] += Field("source_v_x")
	    #result[1] += Field("source_v_x")

        return result

    def bind(self, discr):
        from hedge.mesh import check_bc_coverage
        check_bc_coverage(discr.mesh, [
            self.stressfree_tag,
            self.fixed_tag,
            self.open_tag])

        compiled_op_template = discr.compile(self.op_template())

        def rhs(t, q):
            extra_kwargs = {"q": q}

            if self.source is not None:
                extra_kwargs["source_v_x"] = self.source.volume_interpolant(t, discr)

            return compiled_op_template(**extra_kwargs)

        return rhs

    def max_eigenvalue(self, t, fields=None, discr=None):
	from math import sqrt

        return sqrt(self.C[0,0]/self.rho)

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

