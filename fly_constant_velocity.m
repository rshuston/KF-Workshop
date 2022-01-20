# Constants used here and there

two_pi = 2 * pi;

# Generate truth data

t = [0 : 1 : 100];

n0 = 40;
vn = -1;

e0 = 50;
ve = -1;

n = n0 + vn * t;
e = e0 + ve * t;

r = sqrt(n .* n + e .* e);
theta = atan2(e, n);

rz = r + 0.5 * 2 * (rand(1, length(t)) - 0.5);
thetaz = theta + 0.2 * 2 * (rand(1, length(t)) - 0.5);
if (thetaz >= pi)
    thetaz -= two_pi;
elseif (thetaz < -pi)
    thetaz += two_pi;
endif

nz = rz .* cos(thetaz);
ez = rz .* sin(thetaz);

# Generate filtered data

#[nf, vnf, ef, vef, Pnef] = kf_cv(t, rz, thetaz, 0.1, 0.1, 0.1);
[nf, vnf, ef, vef, Pnef] = ekf_cv(t, rz, thetaz, 0.1, 0.1, 0.1);

# Plot results

rf = sqrt(nf .* nf + ef .* ef);
thetaf = atan2(ef, nf);

figure("name", "Parametric");
plot(ez, nz, ".", ef, nf, "*", e, n, "o");
axis([-50, 50, -50, 50], "square");
set(gca(), "xgrid", "on", "ygrid", "on");

figure("name", "Range");
plot(t, rz, ".", t, rf, "*", t, r, "o");

figure("name", "Angle");
plot(t, thetaz, ".", t, thetaf, "*", t, theta, "o");
