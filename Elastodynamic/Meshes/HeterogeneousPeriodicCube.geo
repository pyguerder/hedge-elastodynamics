// Gmsh project created on Fri May 13 15:55:57 2011
lc = 450;

// Sommet du cube
Point(1) = {1600, -1600, 1600, lc};
Point(2) = {-1600, -1600, 1600, lc};
Point(3) = {-1600, 1600, 1600, lc};
Point(4) = {1600, 1600, 1600, lc};
Point(5) = {1600, -1600, -1600, lc};
Point(6) = {-1600, -1600, -1600, lc};
Point(7) = {-1600, 1600, -1600, lc};
Point(8) = {1600, 1600, -1600, lc};

// Cercles sur la face avant
Point(9) = {-600, 0, 1600, lc};
Point(10) = {0, 0, 1600, lc};
Point(11) = {600, 0, 1600, lc};

// Cercles sur la face arrière
Point(12) = {-600, 0, -1600, lc};
Point(13) = {0, 0, -1600, lc};
Point(14) = {600, 0, -1600, lc};

// Arêtes du cube
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

// Cercle sur la face avant
Circle(13) = {9, 10, 11};
Circle(14) = {11, 10, 9};
Line Loop(39) = {13, 14};
Plane Surface(40) = {39};

// Cercle sur la face arrière
Circle(15) = {12, 13, 14};
Circle(16) = {14, 13, 12};
Line Loop(37) = {15, 16};
Plane Surface(38) = {37};

// Face avant
Line Loop(13) = {1, 2, 3, 4, -13, -14};
Plane Surface(14) = {13};

// Face arrière
Line Loop(15) = {5, 6, 7, 8, -15, -16};
Plane Surface(16) = {15};

// Autres faces
Line Loop(17) = {9, -8, -12, 4};
Plane Surface(18) = {17};
Line Loop(19) = {10, 6, -11, -2};
Plane Surface(20) = {19};
Line Loop(21) = {9, 5, -10, -1};
Plane Surface(22) = {21};
Line Loop(23) = {12, -7, -11, 3};
Plane Surface(24) = {23};

// Limites des deux demi-surfaces de l'inclusion
Line(30) = {9, 12};
Line(32) = {11, 14};

// Demi-surfaces de l'inclusion
Line Loop(33) = {14, 30, -16, -32};
Ruled Surface(34) = {33};
Line Loop(35) = {13, 32, -15, -30};
Ruled Surface(36) = {35};

// Surface et volume du cube moins l'inclusion
Surface Loop(25) = {22, 18, 16, 20, 24, 14, -34, -36};
Volume(43) = {25};

// Surface et volume de l'inclusion
Surface Loop(41) = {34, 40, 36, 38};
Volume(42) = {41};
