function [t, n, e] = fly_wiggle
    tmax = 100;

    t = [0 : 1 : tmax];

    w = 2 * pi / (0.25 * tmax);
    phi0_n = 0;
    phi0_e = pi;
    nw = 2.5;
    ew = 2.5;

    n0 = 40;
    vn = -1;

    e0 = 50;
    ve = -1;

    n = n0 + vn * t + nw * cos(phi0_n + w * t);
    e = e0 + ve * t + ew * cos(phi0_e + w * t);
end
