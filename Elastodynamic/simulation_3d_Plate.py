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
           allow_features=['cuda','mpi'],
           dim=3,
           order=4,
           stfree_tag=TAG_ALL,
           fix_tag=TAG_NONE,
           op_tag=TAG_NONE,
           flux_type="lf",
           max_steps=None,
           output_dir='output',
           #pml=[0, 0, 0, 500, 0, 0],
           pml=None,
           sources=numpy.array([-0.4,-0.3,0.003]),
           source_param={'type': 'Ricker', 'sigma':5e-3, 'fc':20e3, 'td':58e-6, 'amp':1e-2, 'begin':0, 'end':10e-3},
           final_time=40e-3,
           quiet_output=True,
           nonlinearity_type=None,
           mesh_file = 'Meshes/3DPlate.msh',
           periodicity = [None, None, None],
           material_files = ['Materials/aluminium.dat', 'Materials/aluminium.dat'],
           vtu_every=100)

