function fly(scenario)

    # Generate truth data

    if (exist("scenario"))
        switch (scenario)
            case "circle"
                [t, n, e] = fly_circle;
            case "curve"
                [t, n, e] = fly_curve;
            case "line"
                [t, n, e] = fly_line;
            case "square"
                [t, n, e] = fly_square;
            case "wiggle"
                [t, n, e] = fly_wiggle;
            otherwise
                printf("Unknown scenario: %s\n", scenario);
                return;
        endswitch
    else
        [t, n, e] = fly_line;
    endif

    # Generate measurement data

    two_pi = 2 * pi;

    r = sqrt(n .* n + e .* e);
    theta = atan2(e, n);

    rz = r + 0.5 * 2 * (rand(1, length(t)) - 0.5);
    thetaz = theta + 0.2 * 2 * (rand(1, length(t)) - 0.5);
    thetaz = wrapToPi(thetaz);

    nz = rz .* cos(thetaz);
    ez = rz .* sin(thetaz);

    # Generate filtered data

    #[nf, vnf, ef, vef, Pnef] = kf_cv(t, rz, thetaz, 0.1, 0.1, 0.1);
    #[nf, vnf, ef, vef, Pnef] = ekf_cv(t, rz, thetaz, 0.1, 0.1, 0.1);
    [nf, vnf, ef, vef, Pnef] = ukf_cv(t, rz, thetaz, 0.1, 0.1, 0.1);

    # Plot results

    rf = sqrt(nf .* nf + ef .* ef);
    thetaf = atan2(ef, nf);

    figure("name", "Parametric");
    plot(ez, nz, ".", ef, nf, "*", e, n, "o");
    axis([-50, 50, -50, 50], "square");
    set(gca, 'xtick', -50:10:50);
    set(gca, 'ytick', -50:10:50);
    set(gca(), "xgrid", "on", "ygrid", "on");

    figure("name", "Range");
    plot(t, rz, ".", t, rf, "*", t, r, "o");
    set(gca(), "xgrid", "on", "ygrid", "on");

    figure("name", "Angle");
    plot(t, thetaz, ".", t, thetaf, "*", t, theta, "o");
    set(gca(), "xgrid", "on", "ygrid", "on");

endfunction