size_bulk = 0.5;
size_1 = 0.5;
size_2 = 0.5;
size_3 = 0.5;

Point(1) = {-1, -2, 0, size_1};
Point(2) = {-1.3568, -2, 0, size_1};
Point(3) = {-0.6432, -2, 0, size_1};
Point(4) = {-1, -1, 0, size_1};
Point(5) = {-1.3568, -1, 0, size_1};
Point(6) = {-0.6432, -1, 0, size_1};
Point(7) = {-1, 0, 0, size_1};
Point(8) = {-1.3568, 0, 0, size_1};
Point(9) = {-0.6432, 0, 0, size_1};
Point(10) = {-1, 1, 0, size_1};
Point(11) = {-1.3568, 1, 0, size_1};
Point(12) = {-0.6432, 1, 0, size_1};
Point(13) = {-1, 2, 0, size_1};
Point(14) = {-1.3568, 2, 0, size_1};
Point(15) = {-0.6432, 2, 0, size_1};
Point(16) = {0, -2, 0, size_1};
Point(17) = {-0.3568, -2, 0, size_1};
Point(18) = {0.3568, -2, 0, size_1};
Point(19) = {0, -1, 0, size_1};
Point(20) = {-0.3568, -1, 0, size_1};
Point(21) = {0.3568, -1, 0, size_1};
Point(22) = {0, 0, 0, size_3};
Point(23) = {-0.1784, 0, 0, size_3};
Point(24) = {0.1784, 0, 0, size_3};
Point(25) = {0, 1, 0, size_1};
Point(26) = {-0.3568, 1, 0, size_1};
Point(27) = {0.3568, 1, 0, size_1};
Point(28) = {0, 2, 0, size_1};
Point(29) = {-0.3568, 2, 0, size_1};
Point(30) = {0.3568, 2, 0, size_1};
Point(31) = {1, -2, 0, size_1};
Point(32) = {0.6432, -2, 0, size_1};
Point(33) = {1.3568, -2, 0, size_1};
Point(34) = {1, -1, 0, size_1};
Point(35) = {0.6432, -1, 0, size_1};
Point(36) = {1.3568, -1, 0, size_1};
Point(37) = {1, 0, 0, size_1};
Point(38) = {0.6432, 0, 0, size_1};
Point(39) = {1.3568, 0, 0, size_1};
Point(40) = {1, 1, 0, size_1};
Point(41) = {0.6432, 1, 0, size_1};
Point(42) = {1.3568, 1, 0, size_1};
Point(43) = {1, 2, 0, size_1};
Point(44) = {0.6432, 2, 0, size_1};
Point(45) = {1.3568, 2, 0, size_1};
Point(46) = {-6.0, 2.5, 0, size_bulk};
Point(47) = {-6.0, -2.5, 0, size_bulk};
Point(48) = {6.0, 2.5, 0, size_bulk};
Point(49) = {6.0, -2.5, 0, size_bulk};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};
Circle(3) = {5, 4, 6};
Circle(4) = {6, 4, 5};
Circle(5) = {8, 7, 9};
Circle(6) = {9, 7, 8};
Circle(7) = {11, 10, 12};
Circle(8) = {12, 10, 11};
Circle(9) = {14, 13, 15};
Circle(10) = {15, 13, 14};
Circle(11) = {17, 16, 18};
Circle(12) = {18, 16, 17};
Circle(13) = {20, 19, 21};
Circle(14) = {21, 19, 20};
Circle(15) = {23, 22, 24};
Circle(16) = {24, 22, 23};
Circle(17) = {26, 25, 27};
Circle(18) = {27, 25, 26};
Circle(19) = {29, 28, 30};
Circle(20) = {30, 28, 29};
Circle(21) = {32, 31, 33};
Circle(22) = {33, 31, 32};
Circle(23) = {35, 34, 36};
Circle(24) = {36, 34, 35};
Circle(25) = {38, 37, 39};
Circle(26) = {39, 37, 38};
Circle(27) = {41, 40, 42};
Circle(28) = {42, 40, 41};
Circle(29) = {44, 43, 45};
Circle(30) = {45, 43, 44};
Line(31) = {47, 46};
Line(32) = {46, 48};
Line(33) = {48, 49};
Line(34) = {49, 47};
Periodic Line(34) = {-32};
Line Loop(36) = {31, 32, 33, 34, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28, -29, -30};
Plane Surface(37) = {36};
Line Loop(38) = {1, 2};
Plane Surface(39) = {38};
Line Loop(40) = {3, 4};
Plane Surface(41) = {40};
Line Loop(42) = {5, 6};
Plane Surface(43) = {42};
Line Loop(44) = {7, 8};
Plane Surface(45) = {44};
Line Loop(46) = {9, 10};
Plane Surface(47) = {46};
Line Loop(48) = {11, 12};
Plane Surface(49) = {48};
Line Loop(50) = {13, 14};
Plane Surface(51) = {50};
Line Loop(52) = {15, 16};
Plane Surface(53) = {52};
Line Loop(54) = {17, 18};
Plane Surface(55) = {54};
Line Loop(56) = {19, 20};
Plane Surface(57) = {56};
Line Loop(58) = {21, 22};
Plane Surface(59) = {58};
Line Loop(60) = {23, 24};
Plane Surface(61) = {60};
Line Loop(62) = {25, 26};
Plane Surface(63) = {62};
Line Loop(64) = {27, 28};
Plane Surface(65) = {64};
Line Loop(66) = {29, 30};
Plane Surface(67) = {66};

Physical Line("minus_y") = {34};
Physical Line("plus_y") = {32};
Physical Surface("mat1") = {37};
Physical Surface("mat2") = {39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67};

Point(50) = {3, 2.5, 0, size_bulk};
Point(51) = {3, -2.5, 0, size_bulk};
Line(98) = {50, 51};
Physical Line("LineReceivers") = {98};
Line{98} In Surface {37};
