# -*- coding: utf-8 -*-
"""Loads a material into matrixes from a text file"""

from __future__ import division

__copyright__ = "Copyright (C) 2011 Pierre-Yves Guerder <pierre-yves.guerder@centraliens-lille.org>"

__license__ = """
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see U{http://www.gnu.org/licenses/}.
"""

import sys
import numpy
from scipy import triu

class Material:
    # Initialize the material from the given file
    def __init__(self, filename, dtype, print_output):
        
        self.dtype = dtype
        self.print_output = print_output
        self.filename = filename
                
        # open file
        try:
            self.file = open (filename, 'r')
        except IOError, message:
            print >> sys.stderr, 'Failed to open', filename, '-', message
            sys.exit(1);

        # Read the mass density.
        self.rho = self.Read_FloatValue('Density')
        
        # Read the permeability
        self.mu = self.Read_FloatValue('Permeability')
        
        # Read the permittivity
        self.epsilon = self.Read_FloatValue('Permittivity')
        
        # Read the non linearity type
        self.nonlinearity_type = self.Read_TextValue('NonlinearityType')
        
        # Read the linear elastic constants
        self.C, self.elastic_type = self.Read_Matrix('LinearElasticConstants')
        
        # Read the nonlinear elastic constants
        self.Cnl = self.Read_Tensor('NonlinearElasticConstants')
           
        # close file
        self.file.close()
    
    # Read 6 6×6 matrices
    def Read_Tensor(self, name):
        self.file.seek(0, 0)

        # Initialize the Cnl tensor
        Cnl = numpy.zeros((6,6,6), self.dtype)
        
        myline = self.file.readline()
        while ( myline != "" and myline != '$' + name +'\n' ):
            myline = self.file.readline()
    
        try:
            for i in range(6):
                for j in range(6):
                    myline = self.file.readline()
                    if j >= i:
                        myline.strip()
                        elements = myline.split()
                        for k in range(j, 6):
                            Cnl[i][j][k] = float(elements[k])
                            Cnl[j][i][k] = float(elements[k])
                            Cnl[k][j][i] = float(elements[k])
                            Cnl[i][k][j] = float(elements[k])
                            Cnl[k][i][j] = float(elements[k])
                            Cnl[j][k][i] = float(elements[k])
                myline = self.file.readline()
            
        except:
            if self.print_output:
                print self.filename, 'contains no valid', name , 'section.'
            return None, None
        if (myline != '$End' + name + '\n') and self.print_output:
            print self.filename, 'has an invalid', name, 'section.'
        return Cnl
    
    # Read a symmetric 6×6 matrix
    def Read_Matrix(self, name):
        self.file.seek(0, 0)

        # Initialize the C matrix
        C = numpy.zeros((6,6), self.dtype)

        # By default, the elastic type is general
        elastic_type = 'general'
        
        myline = self.file.readline()
        while ( myline != "" and myline != '$' + name +'\n' ):
            myline = self.file.readline()
    
        try:
            for i in range(6):
                myline = self.file.readline()
                myline.strip()
                elements = myline.split()
                for j in range(i, 6):
                    C[i, j] = float(elements[j])
            C = C + triu(C,1).T.astype(self.dtype)
            
            # Check if the elastic type is orthotropic
            temp = True
            for (i, j) in [(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)]:
                temp = temp and C[i, j] == 0
            if(temp):
                elastic_type = 'orthotropic'
            
            # Check if the elastic type is isotropic
            if((C[0, 1] == C[0, 2]) and (C[0, 2] == C[1, 2]) \
                and (C[3, 3] == C[4, 4]) and (C[4, 4] == C[5, 5]) \
                and (C[0, 0] == C[0, 1] + C[3, 3])):
                    elastic_type = 'isotropic'
        except:
            if self.print_output:
                print self.filename, 'contains no valid', name, 'section.'
            return None, None
        if (self.file.readline() != '$End' + name + '\n') and self.print_output:
            print self.filename, 'has an invalid', name, 'section.'
        return C, elastic_type
    
    def Read_TextValue(self, name):
        self.file.seek(0, 0)
        myline = self.file.readline()
        while ( myline != "" and myline != '$' + name +'\n' ):
            myline = self.file.readline()
    
        value = self.file.readline()
        
        if (myline != "" and self.file.readline() != '$End' + name + '\n') and self.print_output:
            print self.filename, 'has an invalid', name, 'section.'
        return value

    def Read_FloatValue(self, name):
        text = self.Read_TextValue(name)
        try:
            value = float(text)
        except ValueError:
            value = 0.
            if self.print_output:
                print self.filename, 'contains no valid', name, 'section.'
            return None
        return value

#myMaterial = Material('Materials/aluminium.dat', numpy.float64)
#print "Rho:", myMaterial.rho
#print "Mu:", myMaterial.mu
#print "Epsilon:", myMaterial.epsilon
#print "C:", myMaterial.C
#print "Elastic type:", myMaterial.elastic_type
#print "NonLinearity type:", myMaterial.nonlinearity_type
#myMaterial = Material('Materials/test.dat', numpy.float64)
#print "Rho:", myMaterial.rho
#print "Mu:", myMaterial.mu
#print "Epsilon:", myMaterial.epsilon
#print "C:", myMaterial.C
#print "Elastic type:", myMaterial.elastic_type
#print "NonLinearity type:", myMaterial.nonlinearity_type
#print "Cnl:", myMaterial.Cnl