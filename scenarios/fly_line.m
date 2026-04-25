function [t, n, e, vn, ve] = fly_line
    tstep = 1;
    t = [0 : tstep : 100];

    n0 = 40;
    vn0 = -1;

    e0 = 50;
    ve0 = -1;

    n = n0 + vn0 * t;
    e = e0 + ve0 * t;
    
    vn = vn0 * ones(size(t));
    ve = ve0 * ones(size(t));
end
