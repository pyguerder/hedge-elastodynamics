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
           allow_features=['mpi'],
           dim=1,
           order=4,
           stfree_tag=TAG_ALL,
           fix_tag=TAG_NONE,
           op_tag=TAG_NONE,
           flux_type="lf",
           max_steps=None,
           output_dir='output',
           pml=None,#[1, 0, 0, 1, 0, 0],
           sources=numpy.array([-7]),
           source_param={'type': 'sinus', 'sigma':0.01, 'fc':722.0, 'td':1/722., 'begin':0, 'end':2/722.},
           final_time=1,
           quiet_output=True,
           nonlinearity_type=None,
           mesh_file = 'Meshes/Resonator1D.msh',
           periodicity = [None],
           material_files = ['Materials/epoxy.dat', 'Materials/aluminium.dat', 'Materials/aluminium.dat'],
           vtu_every=20)

