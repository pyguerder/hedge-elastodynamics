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
           allow_features=['cuda', 'mpi'],
           dim=3,
           order=4,
           stfree_tag=TAG_NONE,
           fix_tag=TAG_ALL,
           op_tag=TAG_NONE,
           flux_type="lf",
           max_steps=None,
           output_dir='output',
           pml=[0, 0, 0, 500, 0, 0],
           sources=numpy.array([1000.0,0.0,0.0]),
           source_param={'type': 'Ricker', 'sigma':10, 'fc':7.25, 'td':0.16, 'begin':0, 'end':2},
           final_time=12,
           quiet_output=True,
           nonlinearity_type=None,
           mesh_file = 'Meshes/HeterogeneousPeriodicCube.msh',
           periodicity = [None, None, None],
           material_files = ['Materials/aluminium.dat', 'Materials/calcite.dat'],
           vtu_every=40)

