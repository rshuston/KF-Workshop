function [t, n, e, vn, ve] = fly_curve
    tstep = 1;
    tmax = 100;

    t = [0 : tstep : tmax];

    phi0 = -pi / 2;
    w = 2 * pi / (4 * tmax);

    nC = -40;
    eC = 40;

    pC = 60;

    n = pC * cos(phi0 + w * t) + nC;
    e = pC * sin(phi0 + w * t) + eC;

    vn = gradient(n, tstep);
    ve = gradient(e, tstep);
end
