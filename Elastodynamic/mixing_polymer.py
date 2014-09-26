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
           order=3,
           stfree_tag=TAG_NONE,
           fix_tag=TAG_ALL,
           op_tag=TAG_NONE,
           flux_type="lf",
           max_steps=25e3,
           output_dir='Mixing2D_polymer_1.3e4_mpi_cq_40Hz',
           pml=[5, 0, 0, 5, 0, 0],
           sources=numpy.array([-19, 0.0]),
           source_param={'type': 'SineBurst', 'sigma':1, 'fc':40, 'td':0.0, 'begin':0, 'end':0.4},
           final_time=1,
           quiet_output=True,
           nonlinearity_type='cubic',
           mesh_file = 'Meshes/Mixing2D.msh',
           periodicity = [None, ('minus_y', 'plus_y')],
           material_files = ['Materials/polymer.dat', 'Materials/polymer.dat'],
           vtu_every=500)

