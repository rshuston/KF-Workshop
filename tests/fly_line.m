function [t, n, e] = fly_line
    t = [0 : 1 : 100];

    n0 = 40;
    vn = -1;

    e0 = 50;
    ve = -1;

    n = n0 + vn * t;
    e = e0 + ve * t;
end
