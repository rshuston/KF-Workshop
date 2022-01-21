function [t, n, e] = fly_square
    ts = [0 : 1 : 24];

    dt = ts(length(ts)) - ts(1);

    topLeft  = [25, -25];
    topRight = [25, 25];
    botRight = [-25, 25];
    botLeft  = [-25, -25];

    dn = topRight(1) - topLeft(1);
    de = topRight(2) - topLeft(2);
    nt = topLeft(1) + dn/dt * ts;
    et = topLeft(2) + de/dt * ts;

    dn = botRight(1) - topRight(1);
    de = botRight(2) - topRight(2);
    nr = topRight(1) + dn/dt * ts;
    er = topRight(2) + de/dt * ts;

    dn = botLeft(1) - botRight(1);
    de = botLeft(2) - botRight(2);
    nb = botRight(1) + dn/dt * ts;
    eb = botRight(2) + de/dt * ts;

    dn = topLeft(1) - botLeft(1);
    de = topLeft(2) - botLeft(2);
    nl = botLeft(1) + dn/dt * ts;
    el = botLeft(2) + de/dt * ts;

    n = [nt, nr, nb, nl];
    e = [et, er, eb, el];
    t = [ts, ts + 25, ts + 50, ts + 75];
endfunction
