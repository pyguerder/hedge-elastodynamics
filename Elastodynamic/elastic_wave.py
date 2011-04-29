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


def main(write_output=True, allow_features='mpi', dim = 2, linear = True,
         stfree_tag=TAG_ALL, fix_tag=TAG_NONE, op_tag=TAG_NONE,
         flux_type_arg="lf", debug=["cuda_no_plan"], dtype = numpy.float64):
    from math import exp
    from libraries.materials import Material
    from hedge.backends import guess_run_context
    rcon = guess_run_context(allow_features)

    if allow_features == 'cuda':
        dtype = numpy.float32
    elif allow_features == 'mpi':
        dtype = numpy.float64

    output_dir = 'output'
    import os
    if not os.access(output_dir, os.F_OK):
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

    if dim == 1:
        if rcon.is_head_rank:
            from hedge.mesh.generator import make_uniform_1d_mesh
            mesh = make_uniform_1d_mesh(-10, 10, 500)
    elif dim == 2:
        if rcon.is_head_rank:
            from hedge.mesh.reader.gmsh import read_gmsh
            mesh_file = 'Meshes/Lamb2Dmod.msh'
            mesh = read_gmsh(mesh_file,
                             force_dimension=2,
                             periodicity=None,
                             allow_internal_boundaries=False,
                             tag_mapper=lambda tag:tag)
    elif dim == 3:
        if rcon.is_head_rank:
            from hedge.mesh.generator import make_ball_mesh
            mesh = make_ball_mesh(max_volume=0.0008)
    else:
        raise RuntimeError, "Bad number of dimensions"


    if rcon.is_head_rank:
        print "%d elements" % len(mesh.elements)
        mesh_data = rcon.distribute_mesh(mesh)
    else:
        mesh_data = rcon.receive_mesh()

    def source_v_x(x, el):
        if dim == 2:
            x = x - numpy.array([1720.0,-2303.28])
        elif dim == 3:
            x = x - numpy.array([1720.0,-2303.28,13.2])
        return exp(-numpy.dot(x, x)*0.01)

    from libraries.functions import TimeRickerWaveletGivenFunction
    from hedge.data import \
            make_tdep_given, \
            TimeIntervalGivenFunction

    discr_init = rcon.make_discretization(mesh_data, order=4)

    # Work out which elements belong to each material
    material_elements = []
    materials2 = []
    for key in materials.keys():
        if key in discr_init.mesh.tag_to_elements.keys():
            material_elements.append(set(discr_init.mesh.tag_to_elements[key]))
        materials2.append(materials[key])

    def mat_val(x, el):
        if len(material_elements) > 0 and el in material_elements[1]:
            return 1
        return 0

    if linear:
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

    discr = rcon.make_discretization(mesh_data, order=4, debug=debug,tune_for=op.op_template())

    from hedge.timestep import LSRK4TimeStepper
    stepper = LSRK4TimeStepper(dtype=dtype)

    from hedge.visualization import VtkVisualizer
    if write_output:
        from os.path import join
        vis = VtkVisualizer(discr, rcon, join(output_dir, "fld"))

    from hedge.tools import join_fields
    fields = join_fields([discr.volume_zeros(dtype=dtype) for _ in range(discr.dimensions)],
            [discr.volume_zeros(dtype=dtype) for _ in range(op.dimF[discr.dimensions])])

    #from hedge.discretization import Filter, ExponentialFilterResponseFunction
    #mode_filter = Filter(discr, ExponentialFilterResponseFunction(min_amplification=0.9, order=4))

    # diagnostics setup -------------------------------------------------------
    from pytools.log import LogManager, \
            add_general_quantities, \
            add_simulation_quantities, \
            add_run_info

    if write_output:
        log_file_name = join(output_dir, "elastic_wave.dat")
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
    rhs = op.bind(discr)
    t=0.0
    try:
        from hedge.timestep import times_and_steps

        step_it = times_and_steps(final_time=2.0,
                                  logmgr=None, #None or logmgr
                                  max_dt_getter=lambda t: op.estimate_timestep(discr, stepper=stepper, t=t, fields=fields))

        for step, t, dt in step_it:
            if step % 50 == 0 and write_output:
                visf = vis.make_file(join(output_dir, "fld-%04d" % step))
                print "%d step" % step

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

if __name__ == "__main__":
    main()
