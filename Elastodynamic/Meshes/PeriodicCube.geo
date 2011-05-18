// Gmsh project created on Fri May 13 15:55:57 2011
lc = 200;

Point(1) = {1600, -1600, 1600, lc};
Point(2) = {-1600, -1600, 1600, lc};
Point(3) = {-1600, 1600, 1600, lc};
Point(4) = {1600, 1600, 1600, lc};
Point(5) = {1600, -1600, -1600, lc};
Point(6) = {-1600, -1600, -1600, lc};
Point(7) = {-1600, 1600, -1600, lc};
Point(8) = {1600, 1600, -1600, lc};


Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};

Line(5) = {7, 8};
Line(6) = {8, 5};
Line(7) = {5, 6};
Line(8) = {6, 7};

Line(9)  = {3, 7};
Line(10) = {4, 8};
Line(11) = {1, 5};
Line(12) = {2, 6};
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};
Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};
Line Loop(17) = {9, -8, -12, 4};
Plane Surface(18) = {17};
Line Loop(19) = {10, 6, -11, -2};
Plane Surface(20) = {19};
Line Loop(21) = {9, 5, -10, -1};
Plane Surface(22) = {21};
Line Loop(23) = {12, -7, -11, 3};
Plane Surface(24) = {23};
Surface Loop(25) = {22, 18, 16, 20, 24, 14};
Volume(26) = {25};
