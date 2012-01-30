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

simulation(write_output=['vtu', 'receivers', 'txt'],
           allow_features=['cuda', 'mpi'],
           dim=2,
           order=5,
           stfree_tag=TAG_NONE,
           fix_tag=TAG_ALL,
           op_tag=TAG_NONE,
           flux_type="lf",
           max_steps=None,
           output_dir='output',
           pml=[400, 400, 0, 400, 400, 0],
           sources=numpy.array([800.0,0.0]),
           source_param={'type': 'Ricker', 'sigma':10, 'fc':7.25, 'td':0.16, 'begin':0, 'end':2},
           final_time=1,
           quiet_output=True,
           nonlinearity_type=None,
           mesh_file = 'Meshes/Rectangle.msh',
           periodicity = [None, None],
           material_files = ['Materials/ex2ddir.dat', 'Materials/ex2ddir.dat', 'Materials/ex2ddir.dat'],
           vtu_every=100)

