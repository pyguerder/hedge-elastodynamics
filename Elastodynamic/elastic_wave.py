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
import numpy.linalg as la
from hedge.mesh import TAG_ALL, TAG_NONE




def main(write_output=True, 
        stfree_tag="stressfree", fix_tag="fixed", op_tag=TAG_NONE, 
        flux_type_arg="lf", debug=["cuda_no_plan"], dtype = numpy.float32):
    from pytools.stopwatch import Job
    from math import sin, cos, pi, exp, sqrt
    from materials import Properties_Mat
    from hedge.backends import guess_run_context
    rcon = guess_run_context(
			     #["mpi"]
			     ["cuda"]
			     )

    dim = 2

    if dim == 1:
        if rcon.is_head_rank:
            from hedge.mesh.generator import make_uniform_1d_mesh
            mesh = make_uniform_1d_mesh(-10, 10, 500)

	rho0, CIJ = Properties_Mat(material="Aluminium",dim=dim,dtype=dtype)	

    elif dim == 2:
        from hedge.mesh.generator import make_rect_mesh
	from hedge.mesh.reader.gmsh import read_gmsh
        if rcon.is_head_rank:
	    mesh = read_gmsh('../../../meshes/gmsh/MeshElastodynamic1.msh', force_dimension=2, periodicity=None,
                              allow_internal_boundaries=False,
                              tag_mapper=lambda tag: tag)

	rho0, CIJ = Properties_Mat(material="Aluminium",dim=dim,dtype=dtype)

    elif dim == 3:
        if rcon.is_head_rank:
            from hedge.mesh.generator import make_ball_mesh
            mesh = make_ball_mesh(max_volume=0.0008)

	rho0, CIJ = Properties_Mat(material="Aluminium",dim=dim,dtype=dtype)
    else:
        raise RuntimeError, "bad number of dimensions"

    if rcon.is_head_rank:
        print "%d elements" % len(mesh.elements)
        mesh_data = rcon.distribute_mesh(mesh)
    else:
        mesh_data = rcon.receive_mesh()

    def source_v_x(x, el):
	#x = x - numpy.array([0.1,0.1])
        return exp(-numpy.dot(x, x)*128)

  # from hedge.models.elastodynamic import ElastoDynamicsOperator
    from elastodynamic import ElastoDynamicsOperator
    from hedge.mesh import TAG_ALL, TAG_NONE
    from hedge.data import \
            make_tdep_given, \
            TimeHarmonicGivenFunction, \
            TimeIntervalGivenFunction

    op = ElastoDynamicsOperator(dimensions=dim, rho=rho0, C=CIJ,
            source=TimeIntervalGivenFunction(
                TimeHarmonicGivenFunction(
                    make_tdep_given(source_v_x), omega=50000),
                0, 0.0001),
            stressfree_tag=stfree_tag,
            fixed_tag=fix_tag,
            open_tag=op_tag,
            flux_type=flux_type_arg,
            )

    discr = rcon.make_discretization(mesh_data, order=4,
		         debug=debug,
			 tune_for=op.op_template()
			            )

    from hedge.timestep import LSRK4TimeStepper
    stepper = LSRK4TimeStepper(dtype=dtype)

    from hedge.visualization import VtkVisualizer
    if write_output:
        vis = VtkVisualizer(discr, rcon, "fld")


    from hedge.tools import join_fields
    fields = join_fields([discr.volume_zeros(dtype=dtype) for i in range(discr.dimensions)],
            [discr.volume_zeros(dtype=dtype) for i in range(op.dimF[discr.dimensions])])

    from hedge.discretization import Filter, ExponentialFilterResponseFunction
    #mode_filter = Filter(discr, ExponentialFilterResponseFunction(min_amplification=0.9, order=4))

    # diagnostics setup -------------------------------------------------------
    from pytools.log import LogManager, \
            add_general_quantities, \
            add_simulation_quantities, \
            add_run_info

    if write_output:
        log_file_name = "elastic_wave.dat"
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

    from hedge.log import Integral, LpNorm
    u_getter = lambda: fields[0]
    logmgr.add_quantity(LpNorm(u_getter, discr, 1, name="l1_u"))
    logmgr.add_quantity(LpNorm(u_getter, discr, name="l2_u"))

    logmgr.add_watches(["step.max", "t_sim.max", "l2_u", "t_step.max"])

    # timestep loop -----------------------------------------------------------
    rhs = op.bind(discr)
    t=0.0
    try:
        from hedge.timestep import times_and_steps

        step_it = times_and_steps(
                final_time=4.0, logmgr=None,
                max_dt_getter=lambda t: op.estimate_timestep(discr,
                    stepper=stepper, t=t, fields=fields))

        for step, t, dt in step_it:
            if step % 10 == 0 and write_output:
                visf = vis.make_file("fld-%04d" % step)
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
