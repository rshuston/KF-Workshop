function [t, n, e] = fly_s_curve
    tmax = 100;

    pC = 25;
    eC = 0;

    t1 = [0 : 1 : tmax/2-1];
    phi0 = 0;
    w = -2 * pi / tmax;
    nC = 25;
    n1 = pC * cos(phi0 + w * t1) + nC;
    e1 = pC * sin(phi0 + w * t1) + eC;

    t2 = [tmax/2 : 1 : tmax];
    w = -w;
    nC = -25;
    n2 = pC * cos(phi0 + w * (t2 - tmax/2)) + nC;
    e2 = pC * sin(phi0 + w * (t2 - tmax/2)) + eC;

    t = [t1, t2];
    n = [n1, n2];
    e = [e1, e2];
end
