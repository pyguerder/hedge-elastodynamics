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

simulation(write_output=True,
           allow_features='mpi',
           dim=2,
           order=4,
           stfree_tag=TAG_NONE,
           fix_tag=TAG_ALL,
           op_tag=TAG_NONE,
           flux_type="lf",
           max_steps=None,
           output_dir='output',
           pml=[400, 400, 0, 400, 400, 0],
           sources=numpy.array([200.0,0.0]),
           final_time=12,
           quiet_output=True,
           nonlinearity_type=None,
           mesh_file = 'Meshes/BiHeterogeneousPeriodicSquare.msh',
           periodicity = [None, None],
           material_files = ['Materials/aluminium.dat', 'Materials/calcite.dat', 'Materials/SiO2.dat'])