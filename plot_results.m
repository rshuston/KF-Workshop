function plot_results(label, t, n, e, r, theta, nz, ez, rz, thetaz, nf, vnf, ef, vef, Pnef)

    rf = sqrt(nf .* nf + ef .* ef);
    thetaf = atan2(ef, nf);

    figure("name", ["Parametric" " - " label]);
    plot(ez, nz, ".", ef, nf, "*", e, n, "o");
    axis([-50, 50, -50, 50], "square");
    set(gca, 'xtick', -50:10:50);
    set(gca, 'ytick', -50:10:50);
    set(gca(), "xgrid", "on", "ygrid", "on");

    figure("name", ["Components" " - " label]);
    
    subplot(2,2,1)
    plot(t, vnf, ".-");
    set(gca(), "xgrid", "on", "ygrid", "on");
    title("vn");

    subplot(2,2,2)
    plot(t, vef, ".-");
    set(gca(), "xgrid", "on", "ygrid", "on");
    title("ve");

    subplot(2,2,3)
    plot(t, rz, ".", t, rf, "*", t, r, "o");
    set(gca(), "xgrid", "on", "ygrid", "on");
    title("Range");

    subplot(2,2,4)
    plot(t, thetaz, ".", t, thetaf, "*", t, theta, "o");
    set(gca(), "xgrid", "on", "ygrid", "on");
    title("Angle");

endfunction