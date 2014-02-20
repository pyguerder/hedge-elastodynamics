# -*- coding: utf-8 -*-

import numpy
from tests import SimulationTestCase
from elastic_wave import main as simulation
from hedge.mesh import TAG_ALL, TAG_NONE

class SimpleTests(SimulationTestCase):
    def test_simulation_2d_linear_plm(self):
        self.clean_folder()
        simulation(write_output=[],
             allow_features=['mpi'],
             dim=2,
             order=4,
             stfree_tag=TAG_NONE,
             fix_tag=TAG_ALL,
             op_tag=TAG_NONE,
             flux_type="lf",
             max_steps=40,
             output_dir=SimulationTestCase.folder,
             pml=[400, 400, 0, 400, 400, 0],
             sources=numpy.array([200.0,0.0]),
             source_param={'type':'sinus', 'sigma':10, 'fc':7.25, 'td':0.16, 'begin':0, 'end':2},
             final_time=12,
             quiet_output=True,
             nonlinearity_type=None,
             mesh_file = 'Meshes/BiHeterogeneousPeriodicSquare.msh',
             periodicity = [None, None],
             material_files = ['Materials/aluminium.dat', 'Materials/calcite.dat', 'Materials/SiO2.dat'],
             vtu_every=20)
        assert True

    def test_simulation_2d_linear_plm_cuda(self):
        self.clean_folder()
        simulation(write_output=[],
             allow_features=['mpi', 'cuda'],
             dim=2,
             order=4,
             stfree_tag=TAG_NONE,
             fix_tag=TAG_ALL,
             op_tag=TAG_NONE,
             flux_type="lf",
             max_steps=40,
             output_dir=SimulationTestCase.folder,
             pml=[400, 400, 0, 400, 400, 0],
             sources=numpy.array([200.0,0.0]),
             source_param={'type':'sinus', 'sigma':10, 'fc':7.25, 'td':0.16, 'begin':0, 'end':2},
             final_time=12,
             quiet_output=True,
             nonlinearity_type=None,
             mesh_file = 'Meshes/BiHeterogeneousPeriodicSquare.msh',
             periodicity = [None, None],
             material_files = ['Materials/aluminium.dat', 'Materials/calcite.dat', 'Materials/SiO2.dat'],
             vtu_every=20)
        assert True

    def test_simulation_3d_nonlinear_nopml(self):
        self.clean_folder()
        simulation(write_output=[],
             allow_features=['mpi'],
             dim=3,
             order=4,
             stfree_tag=TAG_NONE,
             fix_tag=TAG_ALL,
             op_tag=TAG_NONE,
             flux_type="lf",
             max_steps=40,
             output_dir=SimulationTestCase.folder,
             pml=None,
             sources=numpy.array([200.0,0.0,0.0]),
             source_param={'type':'sinus', 'sigma':10, 'fc':7.25, 'td':0.16, 'begin':0, 'end':2},
             final_time=12,
             quiet_output=True,
             nonlinearity_type='classical',
             mesh_file = 'Meshes/HeterogeneousPeriodicCube.msh',
             periodicity = [None, None, None],
             material_files = ['Materials/calcite.dat', 'Materials/calcite.dat', 'Materials/calcite.dat'],
             vtu_every=20)
        assert True

    def test_simulation_2d_linear_plm_mpi(self):
        from pytools.mpi import run_with_mpi_ranks
        run_with_mpi_ranks(__file__, 2, lambda: test_simulation_2d_linear_plm())
