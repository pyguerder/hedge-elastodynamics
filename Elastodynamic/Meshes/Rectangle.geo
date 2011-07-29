lc = 150;

Point(1) = {   0.00,     0.00, 0, lc};
Point(2) = {6000.00,     0.00, 0, lc};
Point(3) = {6000.00, -3000.00, 0, lc};
Point(4) = {   0.00, -3000.00, 0, lc};

// Point sources.
Point(6) = {2000.00,     0.00, 0, lc};

// Receivers.
Point(7) = {3000.00,     0.00, 0, lc};
Point(8) = {4000.00,     0.00, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Point{6, 7, 8} In Surface {5};
Plane Surface(6) = {5};

Physical Point("PointSource") = {6};
Physical Point("PointReceiver") = {6, 7, 8};
Physical Surface("My Surface") = {6};

