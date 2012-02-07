# -*- coding: utf-8 -*-
"""File for defining the parameters of the elastic wave simulation."""

from __future__ import division

__authors__ = ["Olivier Bou Matar <olivier.boumatar@iemn.univ-lille1.fr>",
               "Pierre-Yves Guerder <pierre-yves.guerder@centraliens-lille.org>"]
__copyright__ = "Copyright (C) 2010-2011 the authors"
__license__ = "GNU GPLv3 (or more recent equivalent)"


import numpy
from hedge.mesh import TAG_ALL, TAG_NONE

from elastic_wave import main as simulation

simulation(write_output=['vtu', 'receivers'],
           allow_features=['mpi', 'cuda'],
           dim=2,
           order=4,
           stfree_tag=TAG_NONE,
           fix_tag=TAG_ALL,
           op_tag=TAG_NONE,
           flux_type="lf",
           max_steps=None,
           output_dir='output',
           pml=None,#[0, 0, 0, 0.4, 0, 0],
           sources=numpy.array([-7.0,0.0]),
           source_param={'type': 'sinus', 'sigma':1, 'fc':869.0, 'td':1.1e-3, 'begin':0, 'end':2},
           final_time=1,
           quiet_output=True,
           nonlinearity_type='classical',
           mesh_file = 'Meshes/Resonator5Cells.msh',
           periodicity = [None, ('minus_y', 'plus_y')],
           material_files = ['Materials/epoxy.dat', 'Materials/steel.dat', 'Materials/steel.dat'])

