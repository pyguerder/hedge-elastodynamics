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


def main(write_output=True, allow_features='mpi', dim=2, order=4,
         stfree_tag=TAG_NONE, fix_tag=TAG_ALL, op_tag=TAG_NONE,
         flux_type_arg="lf", debug=[], dtype=numpy.float64,
         max_steps=None, output_dir='output', pml=True,
         override_mesh_sources=False, final_time=12.0, quiet_output=True,
         nonlinearity_type=None):
    """
    Parameters:
    write_output: whether to write (True) visualization files or not (False)
    allow_features: 'mpi' or 'cuda'
    dim: 1, 2 or 3
    order: the order of the method
    stfree_tag: which elements to mark as stress-free boundaries
    fix_tag: which elements to mark as fixed boundaries
    op_tag: which elements to mark as open boundaries
    flux_type: 'lf' (Lax-Freidrich flux) or 'central'
    debug: debug parameters to use in make_discretization()
    dtype: defaults to float64, automatically reduced to float32 for cuda
    max_steps: None (no limit) or maximum number of steps to compute
    output_dir: directory where to write the output
    pml: True or False, to enable or disable the NPML
    override_mesh_sources: if True, ignores the source points of the mesh file
    final_time: number of seconds of simulations to compute
    quiet_output: if True, only the main thread will print information
    nonlinearity_type: None (linear) or 'classical' (non-linear)
    """

    from os import access, makedirs, chdir, F_OK
    from math import exp
    from libraries.materials import Material
    from hedge.backends import guess_run_context
    rcon = guess_run_context(allow_features)
    rcon_init = guess_run_context(allow_features)

    if allow_features == 'cuda':
        dtype = numpy.float32
        debug = ['cuda_no_plan']
    elif allow_features == 'mpi':
        dtype = numpy.float64

    if rcon.is_head_rank and output_dir and not access(output_dir, F_OK):
        makedirs(output_dir)

    if quiet_output:
        print_output = rcon.is_head_rank
    else:
        print_output = True
    if print_output:
        print "Rank", rcon.rank, "will print its output."
    else:
        print "Rank", rcon.rank, "will be silent."

    class Receiver():
        pass

    #
    # Define mesh
    #

    mesh_file = ''

    mesh = None
    if dim == 1:
        from hedge.mesh.generator import make_uniform_1d_mesh
        mesh = make_uniform_1d_mesh(-10, 10, 500)
    elif dim == 2:
        #periodicity= [('minus_x','plus_x'), None]
        periodicity= [None, None]
        mesh_file = 'Meshes/HeterogenPeriodicSquarePML.msh'
        mesh = read_gmsh(mesh_file,
                         force_dimension=2,
                         periodicity=periodicity,
                         allow_internal_boundaries=False,
                         tag_mapper=lambda tag:tag)
    elif dim == 3:
        #periodicity= [('minus_x','plus_x'), ('minus_y','plus_y'), ('minus_z','plus_z')]
        periodicity= [None, None, None]
        mesh_file = 'Meshes/PeriodicCube.msh'
        mesh = read_gmsh(mesh_file,
                         force_dimension=3,
                         periodicity=periodicity,
                         allow_internal_boundaries=False,
                         tag_mapper=lambda tag:tag)
    else:
        raise RuntimeError('Bad number of dimensions')

    if rcon.is_head_rank:
        print "%d elements" % len(mesh.elements)
        mesh_data = rcon.distribute_mesh(mesh)
        mesh_init = rcon_init.distribute_mesh(mesh)
    else:
        mesh_data = rcon.receive_mesh()
        mesh_init = rcon_init.receive_mesh()

    #
    # End of mesh definition
    #

    #
    # Define sources
    #

    source = None
    sources = None
    if mesh_file:
        from libraries.gmsh_reader import GmshReader
        gmsh = GmshReader(mesh_file, dim, print_output)
        sources = gmsh.pointSources
    if sources and not override_mesh_sources:
        #FIXME: Multiple source points are currently unsupported
        source = sources[0]
        if print_output:
            print "Using source from Gmsh file,", source
    else:
        if dim == 2:
            source = numpy.array([800.0,0.0])
        elif dim == 3:
            source = numpy.array([0.0,0.0,0.0])
        if print_output:
            print "Using default source position,", source

    from math import sin, cos, pi

    def source_v_x(x, el):
        x = x - source
        return exp(-numpy.dot(x, x)*0.01)*sin(10*pi/180)

    def source_v_y(y, el):
        y = y - source
        return exp(-numpy.dot(y, y)*0.01)*cos(10*pi/180)

    def source_v_z(z, el):
        z = z - source
        return exp(-numpy.dot(z, z)*0.01)

    from libraries.functions import TimeRickerWaveletGivenFunction
    from hedge.data import make_tdep_given, TimeIntervalGivenFunction

    source_x = TimeIntervalGivenFunction(
                    TimeRickerWaveletGivenFunction(
                        make_tdep_given(source_v_x), fc=7.25, tD = 0.16), 0, 2)
    source_y = TimeIntervalGivenFunction(
                    TimeRickerWaveletGivenFunction(
                        make_tdep_given(source_v_y), fc=7.25, tD = 0.16), 0, 2)
    source_z = TimeIntervalGivenFunction(
                    TimeRickerWaveletGivenFunction(
                        make_tdep_given(source_v_z), fc=7.25, tD = 0.16), 0, 2)

    sources = { 'source_x' : source_x,
                'source_y' : source_y,
                'source_z' : source_z }

    #
    # End of sources definition
    #

    #
    # Define materials and link them with elements
    #

    material1 = Material('Materials/aluminium.dat', dtype, print_output)
    material2 = Material('Materials/calcite.dat', dtype, print_output)

    # Only calcite.dat has a Cnl for the moment!
    if nonlinearity_type is not None:
            material1 = Material('Materials/calcite.dat', dtype, print_output)
            material2 = Material('Materials/calcite.dat', dtype, print_output)

    # Work out which elements belong to each material
    material_elements = []
    materials = []
    speeds = []

    # Default material is material1
    speeds.append((material1.C[0,0]/material1.rho)**0.5)
    materials.append(material1)

    if 'mat1' in mesh_init.tag_to_elements.keys():
        elements_list = [el.id for el in mesh_init.tag_to_elements['mat1']]
        material_elements.append(elements_list)
    if 'mat2' in mesh_init.tag_to_elements.keys():
        elements_list = [el.id for el in mesh_init.tag_to_elements['mat2']]
        material_elements.append(elements_list)
        speeds.append((material2.C[0,0]/material2.rho)**0.5)
        materials.append(material2)
    else:
        # If we have no 'mat2', then the second material is material1
        materials.append(material1)
    speed = max(speeds)

    def mat_val(x, el):
        # Will be used in IfPositive(mat, then, else)
        # 1 will lead to then, 0 to else; default is 0/else
        if len(material_elements) > 1 and el.id in material_elements[1]:
            return 1
        return 0

    #
    # End of materials definition
    #

    #
    # Define the operator
    #

    kwargs = {
              'dimensions': dim,
              'speed': speed,
              'material': make_tdep_given(mat_val),
              'sources': sources,
              'boundaryconditions_tag': \
                    { 'stressfree' : stfree_tag,
                      'fixed' : fix_tag,
                      'open' : op_tag },
              'materials': materials,
              'flux_type': flux_type_arg
              }

    operator = None
    if nonlinearity_type is not None:
        kwargs['nonlinearity_type'] = nonlinearity_type
        if pml:
            from elastodynamic import NLNPMLElastoDynamicsOperator
            operator = 'NLNPMLElastoDynamicsOperator'
        else:
            from elastodynamic import NLElastoDynamicsOperator
            operator = 'NLElastoDynamicsOperator'
    else:
        if pml:
            from elastodynamic import NPMLElastoDynamicsOperator
            operator = "NPMLElastoDynamicsOperator"
        else:
            from elastodynamic import ElastoDynamicsOperator
            operator = "ElastoDynamicsOperator"

    assert operator is not None, "Failed to define operator!"
    op = locals()[operator](**kwargs)

    #
    # End of operator definition
    #

    discr = rcon.make_discretization(mesh_data, order=order, debug=debug,tune_for=op.op_template())

    receivers = None
    point_receivers = []
    i = 0
    if mesh_file:
        gmsh = GmshReader(mesh_file, dim, print_output)
        receivers = gmsh.pointReceivers
    if receivers is not None:
        for receiver in receivers:
            try:
                point_receiver = Receiver()
                point_receiver.evaluator = discr.get_point_evaluator(numpy.array(receiver))
                point_receiver.done_dt = False
                point_receiver.id = i
                point_receiver.coordinates = receiver
            except:
                if print_output:
                    print "Receiver ignored (point not found):", receiver
            else:
                point_receivers.append(point_receiver)
                i += 1
                if print_output:
                    print "Add receiver %d:" % i, receiver

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
        chdir(output_dir)

    if write_output:
        log_file_name = 'elastic_wave.dat'
        for point_receiver in point_receivers:
            point_receiver.pointfile = open("receiver_%s.txt" % point_receiver.id, "wt")
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
        pml_widths = [400, 400, 0, 400, 400, 0]
        coefficients = op.coefficients_from_width(discr, mesh,
                            widths=pml_widths, material=material1, alpha_magnitude=2*pi*0.7)
        rhs = op.bind(discr, coefficients)
    else:
        rhs = op.bind(discr)

    t = 0.0
    max_txt = ''
    try:
        from hedge.timestep import times_and_steps

        step_it = times_and_steps(final_time=final_time,
                                  logmgr=None, #None or logmgr
                                  max_dt_getter=lambda t: op.estimate_timestep(discr, stepper=stepper, t=t, fields=fields))

        for step, t, dt in step_it:
            if max_steps > 0:
                max_txt = ' on ' + format(max_steps)
                if step > max_steps:
                    break

            if step % 20 == 0:
                if print_output:
                    print 'Step ' + format(step) + max_txt

                if write_output:
                    visf = vis.make_file("fld-%04d" % step)
                    variables = [("v", discr.convert_volume(fields[0:dim], "numpy")),
                                 ("F", discr.convert_volume(fields[dim:dim+op.dimF[dim]], "numpy"))]
                    if pml:
                        f2 = ("F2", discr.convert_volume(fields[dim+op.dimF[dim]:dim+op.dimF[dim]+dim*dim*2], "numpy"))
                        variables.append(f2)
                    vis.add_data(visf, variables, time=t, step=step)
                    visf.close()

            if write_output:
                for point_receiver in point_receivers:
                    val = point_receiver.evaluator(fields)
                    buffer = ""
                    if not point_receiver.done_dt:
                        buffer += "# dt: %g s\n" % dt
                        buffer += "# v: %d fields\n" % dim
                        buffer += "# F: %d fields\n" % op.dimF[dim]
                        if pml:
                            buffer += "# F2: %d fields\n" % (dim*dim*2)
                        buffer += "# Coordinates: %s\n" % repr(point_receiver.coordinates)
                        point_receiver.done_dt = True
                    buffer += "%s\n" % str(val)
                    point_receiver.pointfile.write(buffer)

            fields = stepper(fields, t, dt, rhs)
            #fields = mode_filter(fields)

    finally:
        if write_output:
            vis.close()
            for point_receiver in point_receivers:
                point_receiver.pointfile.close()

        logmgr.close()
        discr.close()
        if output_dir:
            chdir('..')

if __name__ == "__main__":
    main()
