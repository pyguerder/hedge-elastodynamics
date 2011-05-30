# -*- coding: utf-8 -*-
"""Function for launching the linear and nonlinear Elastodynamics operators."""

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
from hedge.mesh import TAG_ALL, TAG_NONE
from hedge.mesh.reader.gmsh import read_gmsh

def main(write_output=True, allow_features='mpi', dim = 2, linear = True,
         stfree_tag=TAG_ALL, fix_tag=TAG_NONE, op_tag=TAG_NONE, order = 4,
         flux_type_arg="lf", debug=["cuda_no_plan"], dtype = numpy.float64,
         max_steps = None, output_dir = 'output', pml = False,
         override_mesh_sources = False):
    from math import exp
    from libraries.materials import Material
    from hedge.backends import guess_run_context
    rcon = guess_run_context(allow_features)
    rcon_init = guess_run_context(allow_features)

    if allow_features == 'cuda':
        dtype = numpy.float32
    elif allow_features == 'mpi':
        dtype = numpy.float64

    import os
    if output_dir and not os.access(output_dir, os.F_OK):
        os.makedirs(output_dir)

    mesh_file = ''

    material1 = Material('Materials/aluminium.dat',dtype=dtype)
    material2 = Material('Materials/calcite.dat',dtype=dtype)

    # Only calcite.dat has a Cnl for the moment!
    if linear == False:
            material1 = Material('Materials/calcite.dat',dtype=dtype)
            material2 = Material('Materials/calcite.dat',dtype=dtype)

    materials = { 'mat1': material1,
                  'mat2': material2 }

    speed1 = (material1.C[0,0]/material1.rho)**0.5
    speed2 = (material2.C[0,0]/material2.rho)**0.5
    speed = max(speed1,speed2)

    mesh = None
    if dim == 1:
        if rcon.is_head_rank:
            from hedge.mesh.generator import make_uniform_1d_mesh
            mesh = make_uniform_1d_mesh(-10, 10, 500)
    elif dim == 2:
        #periodicity= [('minus_x','plus_x'), None]
        periodicity= [None, None]
        if rcon.is_head_rank:
            mesh_file = 'Meshes/HeterogeneousPeriodicSquare.msh'
            mesh = read_gmsh(mesh_file,
                             force_dimension=2,
                             periodicity=periodicity,
                             allow_internal_boundaries=False,
                             tag_mapper=lambda tag:tag)
    elif dim == 3:
        #periodicity= [('minus_x','plus_x'), ('minus_y','plus_y'), ('minus_z','plus_z')]
        periodicity= [None, None, None]
        if rcon.is_head_rank:
            mesh_file = 'Meshes/PeriodicCube.msh'
            mesh = read_gmsh(mesh_file,
                             force_dimension=3,
                             periodicity=periodicity,
                             allow_internal_boundaries=False,
                             tag_mapper=lambda tag:tag)
    else:
        raise RuntimeError, "Bad number of dimensions"


    if rcon.is_head_rank:
        print "%d elements" % len(mesh.elements)
        mesh_data = rcon.distribute_mesh(mesh)
        mesh_init = rcon_init.distribute_mesh(mesh)
    else:
        mesh_data = rcon.receive_mesh()
        mesh_init = rcon_init.receive_mesh()

    source = None
    sources = None
    if mesh_file:
        from libraries.gmsh_reader import GmshReader
        gmsh = GmshReader(mesh_file, dim)
        sources = gmsh.pointSources
    if sources and not override_mesh_sources:
        #FIXME: Multiple source points are currently unsupported
        source = sources[0]
        print "Using source from Gmsh file,", source
    else:
        if dim == 2:
            source = numpy.array([800.0,0.0])
        elif dim == 3:
            source = numpy.array([0.0,0.0,0.0])
        print "Using default source position,", source

    def source_v_x(x, el):
        x = x - source
        return exp(-numpy.dot(x, x)*0.01)

    from libraries.functions import TimeRickerWaveletGivenFunction
    from hedge.data import \
            make_tdep_given, \
            TimeIntervalGivenFunction

    # Work out which elements belong to each material
    material_elements = []
    materials2 = []
    for key in materials.keys():
        if key in mesh_init.tag_to_elements.keys():
            elements_list = [el.id for el in mesh_init.tag_to_elements[key]]
            material_elements.append(elements_list)
        materials2.append(materials[key])

    def mat_val(x, el):
        if len(material_elements) > 1 and el.id in material_elements[1]:
                return 1
        return 0

    if linear:
        if pml:
            from elastodynamic import NPMLElastoDynamicsOperator
            op = NPMLElastoDynamicsOperator(dimensions=dim,
                    speed=speed,
                    material = make_tdep_given(mat_val),
                    source=TimeIntervalGivenFunction(
                        TimeRickerWaveletGivenFunction(
                        make_tdep_given(source_v_x), fc=7.25, tD = 0.16),
                        0, 2),
                    boundaryconditions_tag = \
                            { 'stressfree' : stfree_tag,
                              'fixed' : fix_tag,
                              'open' : op_tag },
                    materials = materials2,
                    flux_type=flux_type_arg)
        else:
            from elastodynamic import ElastoDynamicsOperator
            op = ElastoDynamicsOperator(dimensions=dim,
                    speed=speed,
                    material = make_tdep_given(mat_val),
                    source=TimeIntervalGivenFunction(
                        TimeRickerWaveletGivenFunction(
                        make_tdep_given(source_v_x), fc=7.25, tD = 0.16),
                        0, 2),
                    boundaryconditions_tag = \
                            { 'stressfree' : stfree_tag,
                              'fixed' : fix_tag,
                              'open' : op_tag },
                    materials = materials2,
                    flux_type=flux_type_arg)
    else:
        if pml:
            from elastodynamic import NLNPMLElastoDynamicsOperator
            op = NLNPMLElastoDynamicsOperator(dimensions=dim,
                    speed=speed,
                    material = make_tdep_given(mat_val),
                    source=TimeIntervalGivenFunction(
                        TimeRickerWaveletGivenFunction(
                        make_tdep_given(source_v_x), fc=7.25, tD = 0.16),
                        0, 2),
                    nonlinearity_type="classical",
                    boundaryconditions_tag = \
                       { 'stressfree' : stfree_tag,
                         'fixed' : fix_tag,
                         'open' : op_tag },
                    materials = materials2,
                    flux_type=flux_type_arg)
        else:
            from elastodynamic import NLElastoDynamicsOperator
            op = NLElastoDynamicsOperator(dimensions=dim,
                    speed=speed,
                    material = make_tdep_given(mat_val),
                    source=TimeIntervalGivenFunction(
                        TimeRickerWaveletGivenFunction(
                        make_tdep_given(source_v_x), fc=7.25, tD = 0.16),
                        0, 2),
                    nonlinearity_type="classical",
                    boundaryconditions_tag = \
                       { 'stressfree' : stfree_tag,
                         'fixed' : fix_tag,
                         'open' : op_tag },
                    materials = materials2,
                    flux_type=flux_type_arg)


    discr = rcon.make_discretization(mesh_data, order=order, debug=debug,tune_for=op.op_template())

    from hedge.timestep import LSRK4TimeStepper
    stepper = LSRK4TimeStepper(dtype=dtype)

    from hedge.visualization import VtkVisualizer
    if write_output:
        vis = VtkVisualizer(discr, rcon, 'fld')

    from hedge.tools import join_fields
    dim = discr.dimensions
    if pml:
        fields = join_fields([discr.volume_zeros(dtype=dtype) for _ in range(dim)],
                             [discr.volume_zeros(dtype=dtype) for _ in range(op.dimF[dim])],
                             [discr.volume_zeros(dtype=dtype) for _ in range(dim*dim*2)])
    else:
        fields = join_fields([discr.volume_zeros(dtype=dtype) for _ in range(dim)],
                             [discr.volume_zeros(dtype=dtype) for _ in range(op.dimF[dim])])

    #from hedge.discretization import Filter, ExponentialFilterResponseFunction
    #mode_filter = Filter(discr, ExponentialFilterResponseFunction(min_amplification=0.9, order=order))

    # diagnostics setup -------------------------------------------------------
    from pytools.log import LogManager, \
            add_general_quantities, \
            add_simulation_quantities, \
            add_run_info

    if output_dir:
        os.chdir(output_dir)

    if write_output:
        log_file_name = 'elastic_wave.dat'
    else:
        log_file_name = None

    logmgr = LogManager(log_file_name, "w", rcon.communicator)
    add_run_info(logmgr)
    add_general_quantities(logmgr)
    add_simulation_quantities(logmgr)
    discr.add_instrumentation(logmgr)

    from pytools.log import IntervalTimer
    vis_timer = IntervalTimer("t_vis", "Time spent visualizing")
    logmgr.add_quantity(vis_timer)
    stepper.add_instrumentation(logmgr)

    from hedge.log import LpNorm
    u_getter = lambda: fields[0]
    logmgr.add_quantity(LpNorm(u_getter, discr, 1, name="l1_u"))
    logmgr.add_quantity(LpNorm(u_getter, discr, name="l2_u"))

    logmgr.add_watches(["step.max", "t_sim.max", "l2_u", "t_step.max"])

    # timestep loop -----------------------------------------------------------
    if pml:
        # widths: [x_l, y_l, z_l, x_r, y_r, z_r]
        pml_widths = [0, 400, 0, 0, 400, 0]
        coefficients = op.coefficients_from_width(discr, widths=pml_widths, material=material1)
        rhs = op.bind(discr, coefficients)
    else:
        rhs = op.bind(discr)
    
    t=0.0
    max_txt = ''
    try:
        from hedge.timestep import times_and_steps

        step_it = times_and_steps(final_time=2.0,
                                  logmgr=None, #None or logmgr
                                  max_dt_getter=lambda t: op.estimate_timestep(discr, stepper=stepper, t=t, fields=fields))

        for step, t, dt in step_it:
            if max_steps > 0:
                max_txt = ' on ' + format(max_steps)
                if step > max_steps:
                    break
            if step % 20 == 0 and write_output:
                visf = vis.make_file("fld-%04d" % step)
                print 'Step ' + format(step) + max_txt

                vis.add_data(visf,
                        [
                            ("v", discr.convert_volume(fields[0:dim], "numpy")),
                            ("F", discr.convert_volume(fields[dim:], "numpy")),
                        ],
                        time=t,
                        step=step)
                visf.close()

            fields = stepper(fields, t, dt, rhs)
            #fields = mode_filter(fields)

        assert discr.norm(fields) < 1
        assert fields[0].dtype == dtype

    finally:
        if write_output:
            vis.close()

        logmgr.close()
        discr.close()
        if output_dir:
            os.chdir('..')

if __name__ == "__main__":
    main()
