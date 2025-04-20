function fly(scenario)
    
    #
    # NOTE: Coordinate system is NED so that heading is clockwise from north
    #
    
    #
    # Generate truth data
    #
    
    if (exist("scenario"))
        switch (scenario)
            case "circle"
                [t, n, e] = fly_circle;
            case "curve"
                [t, n, e] = fly_curve;
            case "line"
                [t, n, e] = fly_line;
            case "s-curve"
                [t, n, e] = fly_s_curve;
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
    
    #
    # Generate measurement data
    #
    
    two_pi = 2 * pi;
    
    r = sqrt(n .* n + e .* e);
    theta = atan2(e, n); # Heading is clockwise from north
    
    rz = r + 0.1 * 2 * (rand(1, length(t)) - 0.5);
    thetaz = theta + 0.1 * 2 * (rand(1, length(t)) - 0.5);
    thetaz = wrap_minus_plus_pi(thetaz);
    
    ##### Uncomment for noise-free measurements #####
    #rz = r;
    #thetaz = theta;
    #####
    
    nz = rz .* cos(thetaz);
    ez = rz .* sin(thetaz);
    
    #
    # Run and plot Kalman filters
    #
    
    # Linear Kalman filter
    run_filter("KF", @kf_cv_init, @kf_cv_update,
               [10.0, 10.0], [0.1], [0.1, 0.1],
               t, n, e, r, theta, nz, ez, rz, thetaz)
    
    # Extended Kalman filter
    run_filter("EKF", @ekf_cv_init, @ekf_cv_update,
               [10.0, 10.0], [0.1], [0.1, 0.1],
               t, n, e, r, theta, nz, ez, rz, thetaz)
    run_filter("EKF-DC", @ekf_cvdc_init, @ekf_cvdc_update,
               [10.0, 10.0], [0.1], [0.1, 0.1],
               t, n, e, r, theta, nz, ez, rz, thetaz)
    
    # Unscented Kalman filter
#    run_filter("UKF", @ukf_cv_init, @ukf_cv_update,
#               [10.0, 10.0], [0.1], [0.1, 0.1],
#               t, n, e, r, theta, nz, ez, rz, thetaz)
#    run_filter("UKF-DC", @ukf_cvdc_init, @ukf_cvdc_update,
#               [10.0, 10.0], [0.1], [0.1, 0.1],
#               t, n, e, r, theta, nz, ez, rz, thetaz)
#    run_filter("UKF-JU-DC", @ukf_ju_cvdc_init, @ukf_ju_cvdc_update,
#               [10.0, 10.0], [0.1], [0.1, 0.1],
#               t, n, e, r, theta, nz, ez, rz, thetaz)
    
    run_filter("EKF-CT", @ekf_ct_init, @ekf_ct_update,
               [10.0, 10.0, 1.0], [0.1, 0.01], [0.1, 0.1],
                t, n, e, r, theta, nz, ez, rz, thetaz)
    
    run_filter("EKF-CT-DC", @ekf_ctdc_init, @ekf_ctdc_update,
               [10.0, 10.0, 1.0], [0.1, 0.01], [0.1, 0.1],
                t, n, e, r, theta, nz, ez, rz, thetaz)
    
endfunction


function run_filter(label, filter_init, filter_update,
                    init_vars, proc_vars, meas_vars,
                    t, n, e, r, theta, nz, ez, rz, thetaz)
    
    n_kf   = zeros(1, length(t));
    vn_kf  = zeros(1, length(t));
    e_kf   = zeros(1, length(t));
    ve_kf  = zeros(1, length(t));
    Pne_kf = cell(1, length(t));
    
    # k = 1
    s_k = filter_init(t(1), rz(1), thetaz(1), init_vars);
    n_kf(1)  = s_k.x(1);
    vn_kf(1) = s_k.x(2);
    e_kf(1)  = s_k.x(3);
    ve_kf(1) = s_k.x(4);
    Pne_kf{1} = s_k.P;
    
    # k = 2...N
    for k = 2:length(t)
        s_k = filter_update(t(k), s_k, rz(k), thetaz(k), proc_vars, meas_vars);
        n_kf(k)  = s_k.x(1);
        vn_kf(k) = s_k.x(2);
        e_kf(k)  = s_k.x(3);
        ve_kf(k) = s_k.x(4);
        Pne_kf{k} = s_k.P;
    endfor
    
    #s_k
    
    plot_results(label, t, n, e, r, theta, nz, ez, rz, thetaz, n_kf, vn_kf, e_kf, ve_kf, Pne_kf);
    
endfunction


function plot_results(label, t, n, e, r, theta, nz, ez, rz, thetaz, nf, vnf, ef, vef, Pnef)
    
    rf = sqrt(nf .* nf + ef .* ef);
    thetaf = atan2(ef, nf);
    
    figure("name", ["Parametric" " - " label]);
    plot(ez, nz, ".", ef, nf, "*", "markersize", 5, e, n, "o", "markersize", 3);
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
    plot(t, rz, ".", t, rf, "*", "markersize", 5, t, r, "o", "markersize", 3);
    set(gca(), "xgrid", "on", "ygrid", "on");
    title("Range");
    
    subplot(2,2,4)
    plot(t, thetaz, ".", t, thetaf, "*", "markersize", 5, t, theta, "o", "markersize", 3);
    set(gca(), "xgrid", "on", "ygrid", "on");
    title("Angle");
    
endfunction
