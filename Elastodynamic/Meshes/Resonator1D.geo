size_bulk = 0.05;

Point(1) = {-10, 0, 0, size_bulk};
Point(2) = {-4, 0, 0, size_bulk};
Point(3) = {-3, 0, 0, size_bulk};
Point(4) = {-2, 0, 0, size_bulk};
Point(5) = {-1, 0, 0, size_bulk};
Point(6) = {1, 0, 0, size_bulk};
Point(7) = {2, 0, 0, size_bulk};
Point(8) = {3, 0, 0, size_bulk};
Point(9) = {4, 0, 0, size_bulk};
Point(10) = {10, 0, 0, size_bulk};
Point(11) = {-7, 0, 0, size_bulk};
Point(12) = {0, 0, 0, size_bulk};
Point(13) = {0.5, 0, 0, size_bulk};
Point(14) = {7, 0, 0, size_bulk};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};

Physical Line("mat1") = {1, 3, 5, 7, 9};
Physical Line("mat2") = {2, 4, 6, 8};
Physical Point("PointReceiver") = {11, 12, 13, 14};

