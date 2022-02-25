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
    
    #####
    #rz = r;
    #thetaz = theta;
    #####
    
    nz = rz .* cos(thetaz);
    ez = rz .* sin(thetaz);
    
    # Generate filtered data
    
    P0 = [ 1 , 0 , 0 , 0 ;
           0 , 1 , 0 , 0 ;
           0 , 0 , 1 , 0 ;
           0 , 0 , 0 , 1 ];
    
    [n_kf, vn_kf, e_kf, ve_kf, Pne_kf] = kf_cv(t, rz, thetaz, P0, 0.1, 0.1, 0.1);
    [n_ekf, vn_ekf, e_ekf, ve_ekf, Pne_ekf] = ekf_cv(t, rz, thetaz, P0, 0.1, 0.1, 0.1);
    [n_ekf_dc, vn_ekf_dc, e_ekf_dc, ve_ekf_dc, Pne_ekf_dc] = ekf_cvdc(t, rz, thetaz, P0, 0.1, 0.1, 0.1);
    [n_ukf, vn_ukf, e_ukf, ve_ukf, Pne_ukf] = ukf_cv(t, rz, thetaz, P0, 0.1, 0.1, 0.1);
    [n_ukf_dc, vn_ukf_dc, e_ukf_dc, ve_ukf_dc, Pne_ukf_dc] = ukf_cvdc(t, rz, thetaz, P0, 0.1, 0.1, 0.1);
    
    # Plot results
    
    plot_results("KF", t, n, e, r, theta, nz, ez, rz, thetaz, n_kf, vn_kf, e_kf, ve_kf, Pne_kf);
    plot_results("EKF", t, n, e, r, theta, nz, ez, rz, thetaz, n_ekf, vn_ekf, e_ekf, ve_ekf, Pne_ekf);
    plot_results("EKF-DC", t, n, e, r, theta, nz, ez, rz, thetaz, n_ekf_dc, vn_ekf_dc, e_ekf_dc, ve_ekf_dc, Pne_ekf_dc);
    plot_results("UKF", t, n, e, r, theta, nz, ez, rz, thetaz, n_ukf, vn_ukf, e_ukf, ve_ukf, Pne_ukf);
    plot_results("UKF-DC", t, n, e, r, theta, nz, ez, rz, thetaz, n_ukf_dc, vn_ukf_dc, e_ukf_dc, ve_ukf_dc, Pne_ukf_dc);
    
endfunction
