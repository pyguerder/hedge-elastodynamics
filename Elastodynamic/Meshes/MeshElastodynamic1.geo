// Points for the square.
Point(1) = { 1.0,  1.0, 0};
Point(2) = {-1.0,  1.0, 0};
Point(3) = { 1.0, -1.0, 0};
Point(4) = {-1.0, -1.0, 0};

// Points for the circle.
Point(5) = {-0.8,  0.4, 0};
Point(6) = {-0.6,  0.4, 0};
Point(7) = {-0.4,  0.4, 0};

// Lines for the square.
Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};

Circle(5) = {5, 6, 7};
Circle(6) = {7, 6, 5};

Line Loop(7) = {1, 2, 3, 4};
Line Loop(8) = {5, 6};

Plane Surface(9) = {7, 8};

// Boundary conditions.
Physical Line("fixed") = {1, 4, 3, 2};
Physical Line("stressfree") = {6, 5};

Physical Surface("Surf") = {9};

Mesh.SaveElementTagType = 2;
