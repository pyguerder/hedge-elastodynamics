// Gmsh project created on Fri May 13 15:55:57 2011
lc = 100;

Point(1) = {1600, -1600, 0, lc*3};
Point(2) = {-1600, -1600, 0, lc*3};
Point(3) = {-1600, 1600, 0, lc*3};
Point(4) = {1600, 1600, 0, lc*3};
Point(5) = {-1100, 800, 0, lc};
Point(6) = {-800, 800, 0, lc};
Point(7) = {-500, 800, 0, lc};
Point(8) = {800, 0, 0, lc};
Point(9) = {-800, 0, 0, lc};
Point(10) = {-800, 250, 0, lc};
Point(11) = {-800, -250, 0, lc};
Point(12) = {1200, -1200, 0, lc};
Point(13) = {-1200, -1200, 0, lc};
Point(14) = {-1200, 1200, 0, lc};
Point(15) = {1200, 1200, 0, lc};

Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};

Circle(5) = {5, 6, 7};
Circle(6) = {7, 6, 5};

Line(7) = {14, 15};
Line(8) = {15, 12};
Line(9) = {12, 13};
Line(10) = {13, 14};

Line Loop(9) = {7, 8, 9, 10, -5, -6};
Plane Surface(9) = {9};

Line Loop(10) = {5, 6};
Plane Surface(10) = {10};

Line Loop(11) = {1, 2, 3, 4, -7, -8, -9, -10};
Plane Surface(11) = {11};

// periodic boundary conditions - slave on left, master on right
Periodic Line {3} = {-1};
Periodic Line {4} = {-2};

Physical Point("PointSource") = {8};
Physical Point("PointReceiver") = {9, 10, 11};
Physical Line("plus_x") = {2};
Physical Line("minus_x") = {4};
Physical Surface("mat1") = {9};
Physical Surface("mat2") = {10};
Physical Surface("pml") = {11};
