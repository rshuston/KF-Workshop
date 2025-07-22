function [t, n, e] = fly_figure_8
    tmax = 200;

    pC = 25;
    eC = 0;

    t1 = [0 : 1 : tmax/4-1];
    phi0 = 0;
    w = 2 * pi / (tmax / 2);
    nC = 25;
    n1 = pC * cos(phi0 + w * t1) + nC;
    e1 = pC * sin(phi0 + w * t1) + eC;

    t2 = [tmax/4 : 1 : 3*tmax/4-1];
    w = -w;
    nC = -25;
    n2 = pC * cos(phi0 + w * (t2 - tmax/4)) + nC;
    e2 = pC * sin(phi0 + w * (t2 - tmax/4)) + eC;

    t3 = [3*tmax/4 : 1 : tmax];
    phi0 = -pi();
    w = -w;
    nC = 25;
    n3 = pC * cos(phi0 + w * (t3 - 3*tmax/4)) + nC;
    e3 = pC * sin(phi0 + w * (t3 - 3*tmax/4)) + eC;

    t = [t1, t2, t3];
    n = [n1, n2, n3];
    e = [e1, e2, e3];
end
