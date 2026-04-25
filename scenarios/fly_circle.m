function [t, n, e, vn, ve] = fly_circle
    tstep = 1;
    tmax = 100;

    t = [0 : tstep : tmax];

    phi0 = 0;
    w = 2 * pi / tmax;

    nC = 5;
    eC = 5;

    pC = 25;

    n = pC * cos(phi0 + w * t) + nC;
    e = pC * sin(phi0 + w * t) + eC;

    vn = gradient(n, tstep);
    ve = gradient(e, tstep);
end
