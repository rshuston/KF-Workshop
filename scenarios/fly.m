function fly(scenario, add_noise)

    % Default is to fly line
    if nargin < 1
        scenario = 'line';
    end

    % Default is to add noise
    if nargin < 2
        add_noise = 1;
    end

    %
    % NOTE: Coordinate system is NED so that heading is clockwise from north
    %

    %
    % Generate truth data
    %

    switch (scenario)
        case 'circle'
            [t, n, e, vn, ve] = fly_circle;
        case 'curve'
            [t, n, e, vn, ve] = fly_curve;
        case 'figure-8'
            [t, n, e, vn, ve] = fly_figure_8;
        case 'line'
            [t, n, e, vn, ve] = fly_line;
        case 's-curve'
            [t, n, e, vn, ve] = fly_s_curve;
        case 'square'
            [t, n, e, vn, ve] = fly_square;
        case 'wiggle'
            [t, n, e, vn, ve] = fly_wiggle;
        otherwise
            fprintf('Unknown scenario: %s\n', scenario)
            return;
    end

    %
    % Generate measurement data
    %

    r = sqrt(n .* n + e .* e);
    vr = (n .* vn + e .* ve) ./ r;
    vr(isnan(vr)) = 0;
    theta = atan2(e, n); % Heading is clockwise from north

    r_sigma = 0.5;
    vr_sigma = 0.1;
    theta_sigma = deg2rad(4);

    r_var = r_sigma * r_sigma;
    vr_var = vr_sigma * vr_sigma;
    theta_var = theta_sigma * theta_sigma;
    dircos_var = 0.01; % 0.1^2

    if add_noise
        rz = r + r_sigma * 2 * (rand(1, length(t)) - 0.5);
        vrz = vr + vr_sigma * 2 * (rand(1, length(t)) - 0.5);
        thetaz = theta + theta_sigma * 2 * (rand(1, length(t)) - 0.5);
        thetaz = wrap_minus_plus_pi(thetaz);
    else
        rz = r;
        vrz = vr;
        thetaz = theta;
    end

    nz = rz .* cos(thetaz);
    ez = rz .* sin(thetaz);

    %
    % Run and plot filters
    %

    init_var_ne = 100.0;
    init_var_vnve = 1.0;
    init_var_w = 1.0;

    % Linear Kalman filter
    if 1
        run_filter('KF', @kf_cv_init, @kf_cv_update,...
                   [init_var_ne, init_var_vnve], 0.1, [r_var, theta_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Extended Kalman filter
    if 1
        run_filter('EKF', @ekf_cv_init, @ekf_cv_update,...
                   [init_var_ne, init_var_vnve], 0.1, [r_var, theta_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Extended Kalman filter, using direction cosines
    if 0
        run_filter('EKF-DC', @ekf_cvdc_init, @ekf_cvdc_update,...
                   [init_var_ne, init_var_vnve], 0.1, [r_var, dircos_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Extended Kalman filter with range rate
    if 1
        run_filter('EKF2', @ekf2_cv_init, @ekf2_cv_update,...
                   [init_var_ne, init_var_vnve], 0.1, [r_var, theta_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Extended Kalman filter with range rate, using direction cosines
    if 0
        run_filter('EKF2-DC', @ekf2_cvdc_init, @ekf2_cvdc_update,...
                   [init_var_ne, init_var_vnve], 0.1, [r_var, dircos_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Unscented Kalman filter
    if 0
        run_filter('UKF', @ukf_cv_init, @ukf_cv_update,...
                   [init_var_ne, init_var_vnve], 0.1, [r_var, theta_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Unscented Kalman filter, using direction cosines
    if 0
        run_filter('UKF-DC', @ukf_cvdc_init, @ukf_cvdc_update,...
                   [init_var_ne, init_var_vnve], 0.1, [r_var, dircos_var, vr_var],...
                   t, n, e, r, theta,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Unscented Kalman filter, using direction cosines, classic Julier and Uhlmann algorithm
    if 0
        run_filter('UKF-JU-DC', @ukf_ju_cvdc_init, @ukf_ju_cvdc_update,...
                   [init_var_ne, init_var_vnve], 0.1, [r_var, dircos_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Extended Kalman filter, constant turn
    if 0
        run_filter('EKF-CT', @ekf_ct_init, @ekf_ct_update,...
                   [init_var_ne, init_var_vnve, init_var_w], [0.01, 0.01], [r_var, theta_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Extended Kalman filter, constant turn, using direction cosines
    if 0
        run_filter('EKF-CT-DC', @ekf_ctdc_init, @ekf_ctdc_update,...
                   [init_var_ne, init_var_vnve, init_var_w], [0.01, 0.01], [r_var, dircos_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Extended Kalman filter with range rate, constant turn
    if 1
        run_filter('EKF2-CT', @ekf2_ct_init, @ekf2_ct_update,...
                   [init_var_ne, init_var_vnve, init_var_w], [0.01, 0.01], [r_var, theta_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

    % Extended Kalman filter with range rate, constant turn, using direction cosines
    if 0
        run_filter('EKF2-CT-DC', @ekf2_ctdc_init, @ekf2_ctdc_update,...
                   [init_var_ne, init_var_vnve, init_var_w], [0.01, 0.01], [r_var, dircos_var, vr_var],...
                   t, n, e, r, theta, vn, ve, vr,...
                   nz, ez, rz, thetaz, vrz)
    end

end


function run_filter(label, filter_init, filter_update,...
                    init_vars, proc_vars, meas_vars,...
                    t, n, e, r, theta, vn, ve, vr,...
                    nz, ez, rz, thetaz, vrz)

    % k = 1
    s_k = filter_init(t(1), rz(1), thetaz(1), vrz(1), init_vars);

    P_dim = size(s_k.P);

    n_kf   = zeros(1, length(t));
    vn_kf  = zeros(1, length(t));
    e_kf   = zeros(1, length(t));
    ve_kf  = zeros(1, length(t));
    Pne_kf = zeros(P_dim(1), P_dim(2), length(t));

    n_kf(1)  = s_k.x(1);
    vn_kf(1) = s_k.x(2);
    e_kf(1)  = s_k.x(3);
    ve_kf(1) = s_k.x(4);
    Pne_kf(:,:,1) = s_k.P;

    % k = 2...N
    for k = 2:length(t)
        s_k = filter_update(t(k), s_k, rz(k), thetaz(k), vrz(k), proc_vars, meas_vars);
        n_kf(k)  = s_k.x(1);
        vn_kf(k) = s_k.x(2);
        e_kf(k)  = s_k.x(3);
        ve_kf(k) = s_k.x(4);
        Pne_kf(:,:,k) = s_k.P;
    end

    plot_results(label, t, n, e, r, theta, vn, ve, vr,...
                 nz, ez, rz, thetaz, vrz,...
                 n_kf, vn_kf, e_kf, ve_kf, Pne_kf);

end


function plot_results(label, t, n, e, r, theta, vn, ve, vr,...
                      nz, ez, rz, thetaz, vrz,...
                      nf, vnf, ef, vef, Pnef)

    rf = sqrt(nf .* nf + ef .* ef);
    vrf = sqrt(nf .* vnf + ef .* vef);
    thetaf = atan2(ef, nf);

    v = Pnef(1,1,:);
    var_n = v(:).';
    v = Pnef(2,2,:);
    var_vn = v(:).';
    v = Pnef(3,3,:);
    var_e = v(:).';
    v = Pnef(4,4,:);
    var_ve = v(:).';

    figure('name', ['Parametric - ', label]);

    plot(e, n, 'k.', ez, nz, 'r.', ef, nf, 'b.')
    axis([-50, 50, -50, 50], 'square')
    xticks([-50:10:50])
    yticks([-50:10:50])
    legend('true', 'measured', 'filtered', 'Location', 'Best')
    grid on

    figure('name', ['State Components - ' label]);

    subplot(2,2,1)
    plot(t, n, 'k.',t, nz, 'r.',t, nf, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'measured', 'filtered', 'Location', 'Best')
    grid on
    title('n')

    subplot(2,2,2)
    plot(t, e, 'k.',t, ez, 'r.',t, ef, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'measured', 'filtered', 'Location', 'Best')
    grid on
    title('e')

    subplot(2,2,3)
    plot(t, vn, 'k.', t, vnf, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'filtered', 'Location', 'Best')
    grid on
    title('v_n')

    subplot(2,2,4)
    plot(t, ve, 'k.', t, vef, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'filtered', 'Location', 'Best')
    grid on
    title('v_e')

    figure('name', ['State Covariances - ' label]);

    subplot(2,2,1)
    plot(t, var_n, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    grid on
    title('var_n')

    subplot(2,2,2)
    plot(t, var_e, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    grid on
    title('var_e')

    subplot(2,2,3)
    plot(t, var_vn, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    grid on
    title('var_{vn}')

    subplot(2,2,4)
    plot(t, var_ve, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    grid on
    title('var_{ve}')

    figure('name', ['Polar Components - ' label]);

    subplot(3,2,1)
    plot(t, n, 'k.',t, nz, 'r.',t, nf, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'measured', 'filtered', 'Location', 'Best')
    grid on
    title('n')

    subplot(3,2,2)
    plot(t, e, 'k.',t, ez, 'r.',t, ef, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'measured', 'filtered', 'Location', 'Best')
    grid on
    title('e')

    subplot(3,2,3)
    plot(t, vn, 'k.', t, vnf, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'filtered', 'Location', 'Best')
    grid on
    title('v_n')

    subplot(3,2,4)
    plot(t, ve, 'k.', t, vef, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'filtered', 'Location', 'Best')
    grid on
    title('v_e')

    subplot(3,2,5)
    plot(t, r, 'k.', t, rz, 'r.', t, rf, 'b.')
    ylim(marginize_axis_limits(ylim, 1.0));
    legend('true', 'measured', 'filtered', 'Location', 'Best')
    grid on
    title('Range')

    subplot(3,2,6)
    plot(t, rad2deg(theta), 'k.', t, rad2deg(thetaz), 'r.', t, rad2deg(thetaf), 'b.')
    legend('true', 'measured', 'filtered', 'Location', 'Best')
    ylim([-180 180])
    yticks([-180:45:180])
    grid on
    title('Angle')

end
