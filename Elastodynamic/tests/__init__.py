# -*- coding: utf-8 -*-
from unittest import TestCase

from elastic_wave import main as simulation
import os 

class SimulationTestCase(TestCase):
    folder = 'output'

    def clean_folder(self):
        for file in os.listdir(SimulationTestCase.folder):
            file_path = os.path.join(SimulationTestCase.folder, file)
            try:
                os.unlink(file_path)
            except Exception, e:
                print e
