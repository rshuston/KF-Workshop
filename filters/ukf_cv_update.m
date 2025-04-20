# This version of the UKF only approximates the measurement relation, and uses
# different weights for mean (Wm) and covariance (Wc) sums. It also manually
# accounts for angles crossing the [-pi, pi) boundary.

function s_k = ukf_cv_update(t, s_km1, r, theta, proc_vars, meas_vars)
    
    var_p = proc_vars(1);
    var_r = meas_vars(1);
    var_theta = meas_vars(2);
    
    T = t - s_km1.t;
    
    x = s_km1.x;
    P = s_km1.P;
    
    N = s_km1.N;
    M = s_km1.M;
    sqrt_N_plus_lambda = s_km1.sqrt_N_plus_lambda;
    Wm = s_km1.Wm;
    Wc = s_km1.Wc;
    
    two_N_plus_1 = 2 * N + 1;
    
    half_pi = pi / 2;
    
    R = [ var_r , 0         ;
          0     , var_theta ];
    
    # Prediction
    
    Phi = [1 , T , 0 , 0;
           0 , 1 , 0 , 0;
           0 , 0 , 1 , T;
           0 , 0 , 0 , 1];
    
    T2 = T * T;        
    T3d2 = T * T * T / 2;
    T4d4 = T * T * T * T / 4;
    
    Q = [ T4d4 , T3d2 , 0    , 0    ;
          T3d2 , T2   , 0    , 0    ;
          0    , 0    , T4d4 , T3d2 ;
          0    , 0    , T3d2 , T2   ] * var_p;
    
    xp = Phi * x;
    Pp = Phi * P * Phi' + Q;
    
    # Compute sigma points and weights
    
    L = sqrt_N_plus_lambda * chol(Pp, "lower");
    
    X = zeros(N, two_N_plus_1);
    X(:,1) = xp;
    for i = 2:N+1
        X(:,i) = xp + L(:,i-1);
        X(:,i + N) = xp - L(:,i-1);
    endfor
    
    # Compute predicted measurement, saving intermediate results
    theta_ref = atan2(xp(3), xp(1));
    farside = (theta_ref > half_pi || theta_ref < -half_pi);
    
    zp = zeros(M,1);
    h_X = zeros(M, two_N_plus_1);
    for i = 1:two_N_plus_1
        X_i = X(:,i);
        Wm_i = Wm(i);
        h_X_i = [ sqrt(X_i(1) * X_i(1) + X_i(3) * X_i(3)) ;
                  atan2(X_i(3), X_i(1))                   ];
        if (farside)
            h_X_i(2) = wrap_zero_two_pi(h_X_i(2));
        endif
        h_X(:,i) = h_X_i;
        zp += Wm_i * h_X_i;
    endfor
    
    # Normalize angle
    zp(2) = wrap_minus_plus_pi(zp(2));
    
    # Compute covariance matrices
    
    Pxz = zeros(N,M);
    Pzz = zeros(M,M);
    for i = 1:two_N_plus_1
        X_i = X(:,i);
        h_X_i = h_X(:,i);
        Wc_i = Wc(i);
        
        dx = X_i - xp;
        dz = h_X_i - zp;
        dz(2) = wrap_minus_plus_pi(dz(2));
        
        Pxz += Wc_i * (dx * dz');
        Pzz += Wc_i * (dz * dz');
    endfor
    
    # Correction
    
    S = Pzz + R;
    K = Pxz * inv(S);
    
    dr = r - zp(1);
    dtheta = wrap_minus_plus_pi(theta - zp(2));
    dz = [ dr ; dtheta ];
    
    x = xp + K * dz;
    P = Pp - K * S * K';
    
    s_k.t = t;
    s_k.x = x;
    s_k.P = P;
    
    # Carry these forward to the next state update
    s_k.N = N;
    s_k.M = M;
    s_k.sqrt_N_plus_lambda = sqrt_N_plus_lambda;
    s_k.Wm = Wm;
    s_k.Wc = Wc;
    
endfunction
