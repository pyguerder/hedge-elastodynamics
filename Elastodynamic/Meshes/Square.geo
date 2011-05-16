// Gmsh project created on Fri May 13 15:55:57 2011
lc = 200;

Point(1) = {1600, -1600, 0, lc};
Point(2) = {-1600, -1600, 0, lc};
Point(3) = {-1600, 1600, 0, lc};
Point(4) = {1600, 1600, 0, lc};
Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
