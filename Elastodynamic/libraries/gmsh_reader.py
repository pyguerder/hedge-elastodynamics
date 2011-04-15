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

class GmshReader:
    def __init__(self, filename, dim):
        # open file
        try:
            self.file = open (filename, 'r')
        except IOError, message:
            print 'Failed to open file: ', message
            raise RuntimeError("Failed to open file")

        if dim == 1 or dim == 2 or dim == 3:
            self.dim = dim
        else:
            raise RuntimeError, "Bad number of dimensions"

        self.tag_name_map = {}
        self.points_map = {}
        self.nodes_map = {}

        self.Read_PhysicalNames()
        self.Read_Points()
        self.Read_Nodes()

        self.pointReceivers = self.Get_Points('PointReceiver')
        self.pointSources = self.Get_Points('PointSource')

        # close file
        self.file.close()

    def Read_PhysicalNames(self):
        self.file.seek(0, 0)
        myline = self.file.readline()
        while ( myline != "" and myline != '$PhysicalNames\n' ):
            myline = self.file.readline()

        if myline == '':
            print 'This file contains no PhysicalName'
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
                print 'Expected quotes around physical name'

            self.tag_name_map[number, dimension] = name[1:-2]
            name_idx +=1

        if name_count+1 != name_idx:
            print 'Unexpected number of physical names found'

    def Read_Points(self):
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

            if type == 15:
                number_of_tags = int(values[2])
                if number_of_tags == 3:
                    physical_entity = int(values[3])
                    # geometrical_entity = int(values[4])
                    # mesh_partition = int(values[5])
                    node_number_list = int(values[6])
                    self.points_map[number, physical_entity] = node_number_list
                elif number_of_tags == 2:
                    physical_entity = int(values[3])
                    # geometrical_entity = int(values[4])
                    node_number_list = int(values[5])
                    self.points_map[number, physical_entity] = node_number_list
                elif number_of_tags == 1:
                    physical_entity = int(values[3])
                    node_number_list = int(values[4])
                    self.points_map[number, physical_entity] = node_number_list
                elif number_of_tags == 0:
                    node_number_list = int(values[3])
                    self.points_map[number, None] = node_number_list
            name_idx +=1

        if name_count+1 != name_idx:
            print 'Unexpected number of elements found'

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
                print 'Unexpected dimension'

        if name_count+1 != name_idx:
            print 'Unexpected number of nodes found'

    def Get_TagId(self, tag_name):
        if tag_name in self.tag_name_map.values():
            tag_id = [k1 for (k1, k2), v in self.tag_name_map.iteritems() if v == tag_name][0]
            return tag_id
        return None

    def Get_Points(self, tag_name):
        tag_id = self.Get_TagId(tag_name)
        if tag_id:
            elements = [v for (k1, k2), v in self.points_map.iteritems() if k2 == tag_id]
            nodes = [v for k, v in self.nodes_map.iteritems() if k in elements]
            return nodes
        else:
            print "No PointReceiver in this file"

#myGmsh = GmshReader('Meshes/Lamb2D.msh', 3)
#print 'PointReceiver:', myGmsh.pointReceivers
#print 'PointSource:', myGmsh.pointSources