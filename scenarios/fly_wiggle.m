function [t, n, e, vn, ve] = fly_wiggle
    tstep = 1;
    tmax = 100;

    t = [0 : tstep : tmax];

    w = 2 * pi / (0.25 * tmax);
    phi0_n = 0;
    phi0_e = pi;
    nw = 2.5;
    ew = 2.5;

    n0 = 40;
    vn0 = -1;

    e0 = 50;
    ve0 = -1;

    n = n0 + vn0 * t + nw * cos(phi0_n + w * t);
    e = e0 + ve0 * t + ew * cos(phi0_e + w * t);

    vn = gradient(n, tstep);
    ve = gradient(e, tstep);
end
