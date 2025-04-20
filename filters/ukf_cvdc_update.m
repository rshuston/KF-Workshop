# This version of the UKF only approximates the measurement relation, and uses
# different weights for mean (Wm) and covariance (Wc) sums. It uses direction
# cosine angle measurements so that we can easily handle angles crossing the
# [-pi, pi) boundary.

function s_k = ukf_cvdc_update(t, s_km1, r, theta, proc_vars, meas_vars)
    
    var_p = proc_vars(1);
    var_r = meas_vars(1);
    var_dc = meas_vars(2);
    
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
    
    z = [ r      ;
          cos(theta) ;
          sin(theta) ];
    
    # Heuristic covariance, but numerically behaved
    R = [ var_r , 0      , 0      ;
          0     , var_dc , 0      ;
          0     , 0      , var_dc ];
    
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
    
    zp = zeros(M,1);
    h_X = zeros(M, two_N_plus_1);
    for i = 1:two_N_plus_1
        X_i = X(:,i);
        Wm_i = Wm(i);
        rho_i = sqrt(X_i(1) * X_i(1) + X_i(3) * X_i(3));
        h_X_i = [ rho_i          ;
                  X_i(1) / rho_i ;
                  X_i(3) / rho_i ];
        h_X(:,i) = h_X_i;
        zp += Wm_i * h_X_i;
    endfor
    
    # Normalize direction cosines adjusted from weighting
    c = zp(2);
    s = zp(3);
    d = sqrt(c * c + s * s);
    zp(2) = c / d;
    zp(3) = s / d;
    
    # Compute covariance matrices
    
    Pxz = zeros(N,M);
    Pzz = zeros(M,M);
    for i = 1:two_N_plus_1
        X_i = X(:,i);
        h_X_i = h_X(:,i);
        Wc_i = Wc(i);
        
        dx = X_i - xp;
        dz = h_X_i - zp;
        
        Pxz += Wc_i * (dx * dz');
        Pzz += Wc_i * (dz * dz');
    endfor
    
    # Correction
    
    S = Pzz + R;
    K = Pxz * inv(S);
    
    dz = z - zp;
    
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
