# -*- coding: utf-8 -*-
"""Function for launching the linear and nonlinear Elastodynamics operators."""

from __future__ import division

__authors__ = ["Olivier Bou Matar <olivier.boumatar@iemn.univ-lille1.fr>",
               "Pierre-Yves Guerder <pierre-yves.guerder@centraliens-lille.org>"]
__copyright__ = "Copyright (C) 2010-2011 the authors"
__license__ = "GNU GPLv3 (or more recent equivalent)"


import numpy
from hedge.mesh import TAG_ALL, TAG_NONE


def main(write_output=True,
         allow_features='',
         dim=2,
         order=4,
         stfree_tag=TAG_NONE,
         fix_tag=TAG_ALL,
         op_tag=TAG_NONE,
         flux_type="lf",
         max_steps=None,
         output_dir='output',
         pml=None,
         sources=None,
         final_time=12,
         quiet_output=True,
         nonlinearity_type=None,
         mesh_file='',
         periodicity=None,
         material_files=None):
    """
    Parameters:
    @param write_output: whether to write (True) visualization files or not (False)
    @param allow_features: 'mpi' or 'cuda'
    @param dim: 1, 2 or 3
    @param order: the order of the method
    @param stfree_tag: which elements to mark as stress-free boundaries
    @param fix_tag: which elements to mark as fixed boundaries
    @param op_tag: which elements to mark as open boundaries
    @param flux_type: 'lf' (Lax-Freidrich flux) or 'central'
    @param max_steps: None (no limit) or maximum number of steps to compute
    @param output_dir: directory where to write the output
    @param pml: None or NPML widths in this order: [x_l, y_l, z_l, x_r, y_r, z_r]
    @param override_mesh_sources: if True, ignores the source points of the mesh file
    @param final_time: number of seconds of simulations to compute
    @param quiet_output: if True, only the main thread will print information
    @param nonlinearity_type: None (linear) or 'classical' (non-linear)
    @param mesh_file: the file to use as a mesh, or '' in 1D
    @param periodicity: the names of the boundaries to stick together, or None
    @param material_files: array, the material files (.dat) to use
    """

    from os import access, makedirs, chdir, F_OK
    from math import exp, sin, cos, pi
    from libraries.materials import Material
    from hedge.backends import guess_run_context
    from hedge.mesh.reader.gmsh import read_gmsh

    rcon = guess_run_context(allow_features)
    rcon_init = guess_run_context(allow_features)

    debug = []
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

    assert dim in [1, 2, 3], 'Bad number of dimensions'

    # Define mesh ---

    mesh = None
    if mesh_file != '':
        mesh = read_gmsh(mesh_file,
                         force_dimension=dim,
                         periodicity=periodicity,
                         allow_internal_boundaries=False,
                         tag_mapper=lambda tag:tag)
    elif dim == 1:
        from hedge.mesh.generator import make_uniform_1d_mesh
        mesh = make_uniform_1d_mesh(-10, 10, 500)
    else:
        raise Exception('Error: No mesh file specified!')

    if rcon.is_head_rank:
        print "%d elements" % len(mesh.elements)
        mesh_data = rcon.distribute_mesh(mesh)
        mesh_init = rcon_init.distribute_mesh(mesh)
    else:
        mesh_data = rcon.receive_mesh()
        mesh_init = rcon_init.receive_mesh()

    if mesh_file:
        from libraries.gmsh_reader import GmshReader
        gmsh = GmshReader(mesh_file, dim, print_output)

    # End of mesh definition ---
    # Define sources ---

    source = None

    if sources is not None:
        #FIXME: "Multiple source points are currently unsupported"
        source = sources
        if print_output:
            print "Using specified source", source
    else:
        if print_output:
            print "No source specified",
        if mesh_file:
            print "trying to find one in", mesh_file
            sources = gmsh.pointSources
            source = sources[0]
            if sources is not None:
                print "Using source", source, "from", mesh_file
            else:
                print "Error: no source!"
        else:
            print "and no mesh file!"
            raise Exception('Error: Could not find any source!')

    def source_v_x(x, el):
        x = x - source
        return 0

    def source_v_y(y, el):
        y = y - source
        return exp(-numpy.dot(y, y)*0.01)  # cos(10*pi/180)

    def source_v_z(z, el):
        z = z - source
        return 0

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

    # End of sources definition ---
    # Define materials and link them with elements ---

    materials = []
    constants = ['Density', 'LinearElasticConstants']
    if nonlinearity_type is not None:
        constants.append('NonlinearElasticConstants')
    for file in material_files:
        material = Material(file, constants, dtype, print_output)
        if nonlinearity_type is not None:
            # In the nonlinear mode, materials MUST have a nonlinear constants
            assert material.Cnl is not None, "Error: Missing nonlinear constants in " + material_files[0]
        materials.append(material)
    assert len(materials) > 0, "Error: You must define at least 1 material."

    # Work out which elements belong to each material
    material_elements = []
    used_materials = []
    speeds = []

    for num, name in [(0, 'mat1'), (1, 'mat2'), (2, 'mat3')]:
        if len(materials) > num:
            if name in mesh_init.tag_to_elements.keys():
                elements_list = [el.id for el in mesh_init.tag_to_elements[name]]
                material_elements.append(elements_list)
        else:
            num = 0
        speeds.append((materials[num].C[0,0]/materials[num].rho)**0.5)
        used_materials.append(materials[num])
        if print_output:
            print "Using", materials[num].filename, "as", name

    speed = max(speeds)
    if print_output:
        print "Using max speed:", speed, "m/s"

    def mat_val(x, el):
        # Will be used in IfPositive(mat, then, else)
        # 1 will lead to then, 0 to else; default is 0/else
        if len(material_elements) > 2 and el.id in material_elements[2]:
            return 2
        elif len(material_elements) > 1 and el.id in material_elements[1]:
            return 1
        return 0

    # End of materials definition ---
    # Define the elastodynamics operator ---

    kwargs = {
              'dimensions': dim,
              'speed': speed,
              'material': make_tdep_given(mat_val),
              'sources': sources,
              'boundaryconditions_tag': \
                    { 'stressfree' : stfree_tag,
                      'fixed' : fix_tag,
                      'open' : op_tag },
              'materials': used_materials,
              'flux_type': flux_type
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

    # End of elastodynamics operator definition ---
    # Define discretization ---

    discr = rcon.make_discretization(mesh_data, order=order, debug=debug,tune_for=op.op_template())

    # End of discretization definition ---
    # Define receivers ---

    receivers = None
    point_receivers = []
    i = 0
    if mesh_file:
        receivers = gmsh.pointReceivers
    if receivers is not None:
        for receiver in receivers:
            try:
                point_receiver = Receiver()
                point_receiver.evaluator = discr.get_point_evaluator(numpy.array(receiver))
                point_receiver.done_dt = False
                point_receiver.id = i
                point_receiver.coordinates = receiver
                point_receiver.filename = "receiver_%s_%s.txt" % (rcon.rank, point_receiver.id)
            except:
                if not quiet_output:
                    print "Receiver ignored (point not found):", receiver
            else:
                point_receivers.append(point_receiver)
                i += 1
                print "Using", point_receiver.filename, "for receiver", receiver

    # End of receivers definition ---
    # Define timestepping and fields ---

    from hedge.timestep import LSRK4TimeStepper
    stepper = LSRK4TimeStepper(dtype=dtype)

    fields = None
    
    def v():
        if fields is not None:
            return fields[0:dim]
        return [discr.volume_zeros(dtype=dtype) for _ in range(dim)]
    
    def f():
        if fields is not None:
            return fields[dim:dim+op.dimF[dim]]
        return [discr.volume_zeros(dtype=dtype) for _ in range(op.dimF[dim])]
    
    def f2():
        if fields is not None:
            return fields[dim+op.dimF[dim]:dim+op.dimF[dim]+dim*dim*2]
        return [discr.volume_zeros(dtype=dtype) for _ in range(dim*dim*2)]

    fields_list = [v(), f()]
    if pml:
        fields_list.append(f2())

    from hedge.tools import join_fields
    fields = join_fields(*fields_list)

    #from hedge.discretization import Filter, ExponentialFilterResponseFunction
    #mode_filter = Filter(discr, ExponentialFilterResponseFunction(min_amplification=0.9, order=order))

    # End of timestepping and fields definition ---
    # Define visualization ---

    from hedge.visualization import VtkVisualizer
    if write_output:
        vis = VtkVisualizer(discr, rcon, 'fld')

    from pytools.log import LogManager, \
            add_general_quantities, \
            add_simulation_quantities, \
            add_run_info

    if output_dir:
        chdir(output_dir)

    if write_output:
        log_file_name = 'elastic_wave.dat'
        for point_receiver in point_receivers:
            point_receiver.pointfile = open(point_receiver.filename, "wt")
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
    logmgr.add_quantity(LpNorm(u_getter, discr, 1, name="l2_u"))

    logmgr.add_watches(["step.max", "t_sim.max", "l2_u", "t_step.max"])

    # End of visualization definition ---
    # Bind the operator to the discretization ---

    if pml is not None:
        coefficients = op.coefficients_from_width(discr, mesh,
                            widths=pml, material=materials[0], alpha_magnitude=2*pi*0.7)
        rhs = op.bind(discr, coefficients)
    else:
        rhs = op.bind(discr)

    # End of operator binding ---
    # Define the timestep loop ---

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
                    variables = [("v", discr.convert_volume(v(), "numpy")),
                                 ("F", discr.convert_volume(f(), "numpy"))]
                    if pml:
                        f_2 = ("F2", discr.convert_volume(f2(), "numpy"))
                        variables.append(f_2)
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
                        buffer += 't '
                        for i in range(dim):
                            buffer += 'v%s ' % i
                        for i in range(op.dimF[dim]):
                            buffer += "F%s " % i
                        if pml:
                            for i in range(dim*dim*2):
                                buffer += "F''%s " % i
                        point_receiver.done_dt = True
                    buffer += "\n%s " % format(t)
                    for i in range(len(val)):
                        buffer += "%s " % format(val[i])
                    buffer += '\n'
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
            
    # End ---

if __name__ == "__main__":
    main()
