// Gmsh project created on Fri May 13 15:55:57 2011
lc = 200;

Point(1) = {1600, -1600, 0, lc};
Point(2) = {-1600, -1600, 0, lc};
Point(3) = {-1600, 1600, 0, lc};
Point(4) = {1600, 1600, 0, lc};
Point(5) = {-1300, 800, 0, lc};
Point(6) = {-800, 800, 0, lc};
Point(7) = {-300, 800, 0, lc};
Point(8) = {800, 0, 0, lc};
Point(9) = {-800, 0, 0, lc};

Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};

Circle(5) = {5, 6, 7};
Circle(6) = {7, 6, 5};

Line Loop(9) = {1, 2, 3, 4, -5, -6};
Plane Surface(9) = {9};

Line Loop(10) = {5, 6};
Plane Surface(10) = {10};

// periodic boundary conditions - slave on left, master on right
Periodic Line {3} = {-1};
Periodic Line {4} = {-2};

Physical Point("PointSource") = {8};
Physical Point("PointReceiver") = {9};
Physical Line("plus_x") = {2};
Physical Line("minus_x") = {4};
Physical Surface("mat1") = {9};
Physical Surface("mat2") = {10};
