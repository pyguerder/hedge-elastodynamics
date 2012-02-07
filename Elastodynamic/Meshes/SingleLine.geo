//SingleLine, created with crystalpy
size_bulk = 0.5;
size_1 = 0.5;
size_2 = 0.5;
size_3 = 0.5;

Point(1) = {-3, 0, 0, size_1};
Point(2) = {-3.3568, 0, 0, size_1};
Point(3) = {-2.6432, 0, 0, size_1};
Point(4) = {-2, 0, 0, size_3};
Point(5) = {-2.1784, 0, 0, size_3};
Point(6) = {-1.8216, 0, 0, size_3};
Point(7) = {-1, 0, 0, size_2};
Point(8) = {-1.1262, 0, 0, size_2};
Point(9) = {-0.8738, 0, 0, size_2};
Point(10) = {-6.5, 0.5, 0, size_bulk};
Point(11) = {-6.5, -0.5, 0, size_bulk};
Point(12) = {6.5, 0.5, 0, size_bulk};
Point(13) = {6.5, -0.5, 0, size_bulk};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};
Circle(3) = {5, 4, 6};
Circle(4) = {6, 4, 5};
Circle(5) = {8, 7, 9};
Circle(6) = {9, 7, 8};
Line(7) = {11, 10};
Line(8) = {10, 12};
Line(9) = {12, 13};
Line(10) = {13, 11};
Periodic Line(10) = {-8};
Line Loop(12) = {7, 8, 9, 10, -1, -2, -3, -4, -5, -6};
Plane Surface(13) = {12};
Line Loop(14) = {1, 2};
Plane Surface(15) = {14};
Line Loop(16) = {3, 4};
Plane Surface(17) = {16};
Line Loop(18) = {5, 6};
Plane Surface(19) = {18};

Physical Line("minus_y") = {10};
Physical Line("plus_y") = {8};
Physical Surface("mat1") = {13};
Physical Surface("mat2") = {15, 17, 19};
