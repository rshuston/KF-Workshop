# This version of the EKF manually accounts for angles crossing the
# [-pi, pi) boundary.

function s_k = ekf_cv_update(t, s_km1, r, theta, proc_vars, meas_vars)
    
    var_p = proc_vars(1);
    var_r = meas_vars(1);
    var_theta = meas_vars(2);
    
    T = t - s_km1.t;
    
    x = s_km1.x;
    P = s_km1.P;
    
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
    
    R = [ var_r , 0         ;
          0     , var_theta ];
    
    xp = Phi * x;
    Pp = Phi * P * Phi' + Q;
    
    n2e2 = xp(1) * xp(1) + xp(3) * xp(3);
    
    rp = sqrt(n2e2);
    thetap = atan2(xp(3), xp(1));
    
    dr = r - rp;
    dtheta = wrap_minus_plus_pi(theta - thetap);
    dz = [ dr ; dtheta ];
    
    H = [ xp(1) / sqrt(n2e2) , 0 , xp(3) / sqrt(n2e2) , 0 ;
          - xp(3) / n2e2 ,     0 , xp(1) / n2e2       , 0 ];
    
    S = H * Pp * H' + R;
    K = Pp * H' * inv(S);
    
    x = xp + K * dz;
    P = Pp - K * S * K';
    
    s_k.t = t;
    s_k.x = x;
    s_k.P = P;
    
endfunction
