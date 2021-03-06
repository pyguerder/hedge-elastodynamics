# -*- coding: utf8 -*-
"""Loads the physical points of a Gmsh file using the Gmsh reader"""

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

class GmshReader:
    def __init__(self, filename, dim, print_output):
        self.filename = filename

        # open file
        try:
            self.file = open (filename, 'r')
        except IOError, message:
            print >> sys.stderr, 'Failed to open file: ', message
            raise RuntimeError("Failed to open file")

        self.print_output = print_output

        if dim == 1 or dim == 2 or dim == 3:
            self.dim = dim
        else:
            raise RuntimeError, "Bad number of dimensions"

        self.tag_name_map = {}
        self.points_map = {}
        self.lines_map = {}
        self.surfaces_map = {}
        self.nodes_map = {}

        self.Read_PhysicalNames()
        self.Read_Elements()
        self.Read_Nodes()

        self.pointReceivers = self.Get_Points('PointReceiver')
        self.pointSources = self.Get_Points('PointSource')
        self.internalBoundaries = self.Get_Lines('MyLine')
        self.pointReceivers.extend(self.Get_Points_of_Lines('LineReceivers'))

        if dim == 1:
            self.materials = \
                {'mat1' : self.Get_Lines('mat1'), 
                 'mat2' : self.Get_Lines('mat2'),
                 'mat3' : self.Get_Lines('mat3')}
        else:
            self.materials = \
                {'mat1' : self.Get_Surfaces('mat1'), 
                 'mat2' : self.Get_Surfaces('mat2'),
                 'mat3' : self.Get_Surfaces('mat3')}

        # close file
        self.file.close()

    def Read_PhysicalNames(self):
        self.file.seek(0, 0)
        myline = self.file.readline()
        while ( myline != "" and myline != '$PhysicalNames\n' ):
            myline = self.file.readline()

        if myline == '':
            if self.print_output:
                print self.filename, 'contains no PhysicalName'
            return

        myline = self.file.readline()
        name_count = int(myline)
        name_idx = 1

        while True:
            myline = self.file.readline()
            if myline == '$EndPhysicalNames\n':
                    break

            dimension, number, name = myline.split(' ', 2)

            dimension = int(dimension)
            number = int(number)

            if not name[0] == '"' or not name[-2] == '"':
                print >> sys.stderr, 'Expected quotes around physical name'

            self.tag_name_map[number, dimension] = name[1:-2]
            name_idx +=1

        if name_count+1 != name_idx:
            print >> sys.stderr, 'Unexpected number of physical names found'

    def Read_Elements(self):
        # 
        # Description of the Elements section:
        # Type 1 (line): 2 tags => 7 values, 3 tags => 8 values:
        # [id] 1 [num tags] [tag1] [tag2] ... [node1] [node2]
        # Type 2 (surface): 2 tags => 8 values, 3 tags => 9 values:
        # [id] 2 [num tags] [tag1] [tag2] ... [node1] [node2] [node3]
        # Type 15 (physical point): 2 tags => 6 values, 3 tags => 7 values:
        # [id] 15 [num tags] [tag1] [tag2] ... [node]
        # etc.
        #
        self.file.seek(0, 0)
        myline = self.file.readline()
        while ( myline != "" and myline != '$Elements\n' ):
            myline = self.file.readline()

        myline = self.file.readline()
        name_count = int(myline)
        name_idx = 1

        while True:
            myline = self.file.readline()
            if myline == '$EndElements\n':
                break

            values = myline.split()

            number = int(values[0])
            type = int(values[1])
            number_of_tags = int(values[2])

            if type == 15:  # 15 means Point
                physical_point = None
                if number_of_tags != 0:
                    physical_point = int(values[3])
                node = int(values[number_of_tags + 3])
                self.points_map[number, physical_point] = node
            elif type == 2:  # 2 means Surface
                physical_surface = None
                if number_of_tags != 0:
                    physical_surface = int(values[3])
                node1 = int(values[number_of_tags + 3])
                node2 = int(values[number_of_tags + 4])
                node3 = int(values[number_of_tags + 5])
                self.surfaces_map[number, physical_surface] = (node1, node2, node3)
            elif type == 1:  # 1 means Line
                physical_line = None
                if number_of_tags != 0:
                    physical_line = int(values[3])
                node1 = int(values[number_of_tags + 3])
                node2 = int(values[number_of_tags + 4])
                self.lines_map[number, physical_line] = (node1, node2)
            name_idx +=1

        if name_count+1 != name_idx:
            print >> sys.stderr, 'Unexpected number of elements found'

    def Read_Nodes(self):
        self.file.seek(0, 0)
        myline = self.file.readline()
        while ( myline != "" and myline != '$Nodes\n' ):
            myline = self.file.readline()

        myline = self.file.readline()
        name_count = int(myline)
        name_idx = 1

        while True:
            myline = self.file.readline()
            if myline == '$EndNodes\n':
                break

            values = myline.split()

            try:
                number = int(values[0])
                if self.dim == 3:
                    coord1 = float(values[1])
                    coord2 = float(values[2])
                    coord3 = float(values[3])
                    self.nodes_map[number] = (coord1, coord2, coord3)
                elif self.dim == 2:
                    coord1 = float(values[1])
                    coord2 = float(values[2])
                    self.nodes_map[number] = (coord1, coord2)
                elif self.dim == 1:
                    coord1 = float(values[1])
                    self.nodes_map[number] = (coord1)
                name_idx +=1
            except:
                print >> sys.stderr, 'Unexpected dimension'

        if name_count+1 != name_idx:
            print >> sys.stderr, 'Unexpected number of nodes found'

    def Get_TagId(self, tag_name):
        if tag_name in self.tag_name_map.values():
            tag_id = [k1 for (k1, _), v in self.tag_name_map.iteritems() if v == tag_name][0]
            return tag_id
        return None

    def Get_Points(self, tag_name):
        tag_id = self.Get_TagId(tag_name)
        if tag_id:
            elements = [v for (_, k2), v in self.points_map.iteritems() if k2 == tag_id]
            nodes = [v for k, v in self.nodes_map.iteritems() if k in elements]
            return nodes
        else:
            if self.print_output:
                print self.filename, 'contains no', tag_name
            return []

    def Get_Lines(self, tag_name):
        tag_id = self.Get_TagId(tag_name)
        if tag_id:
            elements = [k1 for (k1, k2), _ in self.lines_map.iteritems() if k2 == tag_id]
            return elements
        else:
            if self.print_output:
                print self.filename, 'contains no', tag_name
            return []

    def Get_Points_of_Lines(self, tag_name):
        tag_id = self.Get_TagId(tag_name)
        if tag_id:
            lines = [k1 for (k1, k2), _ in self.lines_map.iteritems() if k2 == tag_id]
            nodes = []
            for (k1, k2), v in self.lines_map.iteritems():
                if k1 in lines:
                    nodes.extend(v)
            coords = [v for k, v in self.nodes_map.iteritems() if k in nodes]
            return coords
        else:
            if self.print_output:
                print self.filename, 'contains no', tag_name
            return []

    def Get_Surfaces(self, tag_name):
        tag_id = self.Get_TagId(tag_name)
        if tag_id:
            elements = [k1 for (k1, k2), _ in self.surfaces_map.iteritems() if k2 == tag_id]
            return elements
        else:
            if self.print_output:
                print self.filename, 'contains no', tag_name
            return []

#myGmsh = GmshReader('../Meshes/SquareIntBound.msh', 2, True)
#print 'PointReceiver:', myGmsh.pointReceivers
#print 'PointSource:', myGmsh.pointSources
#print 'Materials:', myGmsh.materials
#print 'Line:', myGmsh.internalBoundaries
