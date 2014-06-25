//Mixing2D, created with crystalpy
size_bulk = 0.8;
size_1 = 0.0818535277187;

Point(1) = {-15, 0, 0, 0.8};
Point(2) = {0, 0, 0, 0.8};
Point(3) = {12, 0, 0, 0.8};
Point(4) = {23, 0, 0, 0.8};
Point(5) = {-23, 0.4, 0, 0.8};
Point(6) = {-23, -0.4, 0, 0.8};
Point(7) = {-0.0409267638594, 0.0409267638594, 0, size_1};
Point(8) = {-0.0409267638594, -0.0409267638594, 0, size_1};
Point(9) = {0.0409267638594, 0.0409267638594, 0, size_1};
Point(10) = {0.0409267638594, -0.0409267638594, 0, size_1};
Point(11) = {-25.0, 0.5, 0, size_bulk};
Point(12) = {-25.0, -0.5, 0, size_bulk};
Point(13) = {25.0, 0.5, 0, size_bulk};
Point(14) = {25.0, -0.5, 0, size_bulk};

Line(1) = {5, 6};
Line(2) = {8, 7};
Line(3) = {7, 9};
Line(4) = {9, 10};
Line(5) = {10, 8};
Line(6) = {12, 11};
Line(7) = {11, 13};
Line(8) = {13, 14};
Line(9) = {14, 12};
Periodic Line(9) = {-7};
Line Loop(11) = {6, 7, 8, 9, -2, -3, -4, -5};
Plane Surface(12) = {11};
Line Loop(13) = {2, 3, 4, 5};
Plane Surface(14) = {13};

Physical Point("PointReceiver") = {1, 2, 3, 4};
Physical Line("LineSource") = {1};
Physical Line("minus_y") = {9};
Physical Line("plus_y") = {7};
Physical Surface("mat1") = {12};
Line{1} In Surface {12};
Physical Surface("mat2") = {14};
