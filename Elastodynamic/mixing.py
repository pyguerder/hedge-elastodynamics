# -*- coding: utf-8 -*-
"""File for defining the parameters of the elastic wave simulation."""

from __future__ import division

__authors__ = ["Olivier Bou Matar <olivier.boumatar@iemn.univ-lille1.fr>",
               "Pierre-Yves Guerder <pierre-yves.guerder@centraliens-lille.org>"]
__copyright__ = "Copyright (C) 2010-2011 the authors"
__license__ = "GNU GPLv3 (or more recent equivalent)"


from hedge.mesh import TAG_ALL, TAG_NONE

from elastic_wave import main as simulation

import numpy

simulation(write_output=['vtu', 'receivers'],
           allow_features=['mpi'],
           dim=2,
           order=4,
           stfree_tag=TAG_NONE,
           fix_tag=TAG_ALL,
           op_tag=TAG_NONE,
           flux_type="lf",
           max_steps=75e3,
           output_dir='output',
           pml=[0.4, 0, 0, 0.4, 0, 0],
           sources=numpy.array([-23.0, 0.0]),
           source_param={'type': 'SineBurst', 'sigma':1, 'fc':200, 'td':0.0, 'begin':0, 'end':0.25},
           final_time=12,
           quiet_output=True,
           nonlinearity_type='cubic',
           mesh_file = 'Meshes/Mixing2D.msh',
           periodicity = [None, ('minus_y', 'plus_y')],
           material_files = ['Materials/PDMS.dat', 'Materials/polymer.dat'],
           vtu_every=100)

