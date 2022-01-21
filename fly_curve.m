function [t, n, e] = fly_curve
    tmax = 100;

    t = [0 : 1 : tmax];

    phi0 = -pi / 2;
    w = 2 * pi / (4 * tmax);

    nC = -40;
    eC = 40;

    pC = 60;

    n = pC * cos(phi0 + w * t) + nC;
    e = pC * sin(phi0 + w * t) + eC;
endfunction
