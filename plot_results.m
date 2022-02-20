function plot_results(label, t, n, e, r, theta, nz, ez, rz, thetaz, nf, vnf, ef, vef, Pnef)

    rf = sqrt(nf .* nf + ef .* ef);
    thetaf = atan2(ef, nf);

    figure("name", ["Parametric" " - " label]);
    plot(ez, nz, ".", ef, nf, "*", e, n, "o");
    axis([-50, 50, -50, 50], "square");
    set(gca, 'xtick', -50:10:50);
    set(gca, 'ytick', -50:10:50);
    set(gca(), "xgrid", "on", "ygrid", "on");

    figure("name", ["vn" " - " label]);
    plot(t, vnf, ".-");
    set(gca(), "xgrid", "on", "ygrid", "on");

    figure("name", ["ve" " - " label]);
    plot(t, vef, ".-");
    set(gca(), "xgrid", "on", "ygrid", "on");

    figure("name", ["Range" " - " label]);
    plot(t, rz, ".", t, rf, "*", t, r, "o");
    set(gca(), "xgrid", "on", "ygrid", "on");

    figure("name", ["Angle" " - " label]);
    plot(t, thetaz, ".", t, thetaf, "*", t, theta, "o");
    set(gca(), "xgrid", "on", "ygrid", "on");

endfunction