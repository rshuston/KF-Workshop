t = [0 : 1 : 1];

n0 = 1;
vn = 1;

e0 = 1;
ve = 0;

n = n0 + vn * t;
e = e0 + ve * t;

r = sqrt(n .* n + e .* e);
theta = atan2(e, n);

nz = r .* cos(theta);
ez = r .* sin(theta);

n
e
r
theta
nz
ez

[nf, vnf, ef, vef, Pnef] = ekf_cv(t, r, theta, 1, 1, 1);

nf
ef
Pnef
