lc_particle = 0.03;
lc_outer = 0.08;
a = 0.2;
R = 1.0;

Point(1) = {0.0, 0.0, 0.0, lc_particle};

Point(2) = { a, 0.0, 0.0, lc_particle};
Point(3) = {0.0,  a, 0.0, lc_particle};
Point(4) = {-a, 0.0, 0.0, lc_particle};
Point(5) = {0.0, -a, 0.0, lc_particle};

Point(6) = { R, 0.0, 0.0, lc_outer};
Point(7) = {0.0,  R, 0.0, lc_outer};
Point(8) = {-R, 0.0, 0.0, lc_outer};
Point(9) = {0.0, -R, 0.0, lc_outer};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {2, 1};

Physical Surface("particle", 1) = {1};
Physical Surface("background", 2) = {2};
Physical Curve("particle_boundary", 11) = {1, 2, 3, 4};
Physical Curve("outer_boundary", 12) = {5, 6, 7, 8};
