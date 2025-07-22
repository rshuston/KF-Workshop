function [t, n, e] = fly_circle
    tmax = 100;

    t = [0 : 1 : tmax];

    phi0 = 0;
    w = 2 * pi / tmax;

    nC = 5;
    eC = 5;

    pC = 25;

    n = pC * cos(phi0 + w * t) + nC;
    e = pC * sin(phi0 + w * t) + eC;
end
