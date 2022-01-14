# Generate truth data

t = [0 : 1 : 100];

n0 = 40;
vn = -1;

e0 = 50;
ve = -1;

n = n0 + vn * t;
e = e0 + ve * t;

# Generate measurement data

#r = sqrt(n .* n + e .* e);
#theta = atan2(e, n);

r = sqrt(n .* n + e .* e) + 0.5 * 2 * (rand(1, length(t)) - 0.5);
theta = atan2(e, n) + 0.2 * 2 * (rand(1, length(t)) - 0.5);

nz = r .* cos(theta);
ez = r .* sin(theta);

# Generate filtered data

#[nf, vnf, ef, vef, Pnef] = kf_cv(t, r, theta, 0.1, 0.1, 0.1);
[nf, vnf, ef, vef, Pnef] = ekf_cv(t, r, theta, 0.1, 0.1, 0.1);

# Plot results

rf = sqrt(nf .* nf + ef .* ef);
thetaf = atan2(ef, nf);

figure("name", "Parametric");
plot(e, n, ".", ez, nz, "o", ef, nf, "*");
axis([-50, 50, -50, 50], "square");
set(gca(), "xgrid", "on", "ygrid", "on");

figure("name", "Range");
plot(t, r, ".", t, rf, "*");

figure("name", "Angle");
plot(t, theta, ".", t, thetaf, "*");
