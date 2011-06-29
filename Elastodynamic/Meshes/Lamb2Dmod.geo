lc = 70;

Point(1) = {   0.00,     0.00, 0, lc};
Point(2) = {4000.00,     0.00, 0, lc};
Point(3) = {4000.00, -2705.31, 0, lc};
Point(4) = {   0.00, -2000.00, 0, lc};

// Point sources.
Point(6) = {1720.00, -2303.28, 0, lc};

// Receivers.
Point(7) = {2557.10, -2450.885, 0, lc};
Point(8) = {2901.80, -2511.665, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Physical Point("PointSource") = {6};
Physical Point("PointReceiver") = {7, 8};
Physical Surface("My Surface") = {6};
