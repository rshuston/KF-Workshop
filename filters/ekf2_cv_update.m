% This version of the EKF manually accounts for angles crossing the
% [-pi, pi) boundary.

function s_k = ekf2_cv_update(t, s_km1, r, theta, vr, proc_vars, meas_vars)
    
    var_p = proc_vars(1);
    var_r = meas_vars(1);
    var_theta = meas_vars(2);
    var_vr = meas_vars(3);
    
    T = t - s_km1.t;
    
    x = s_km1.x;
    P = s_km1.P;
    
    % Prediction
    
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
    
    R = [ var_r , 0      , 0         ;
          0     , var_vr , 0         ;
          0     , 0      , var_theta ];
    
    xp = Phi * x;
    % Pp = Phi * P * Phi' + Q;
    Pp = quad_mult(Phi, P) + Q;
    
    rp2 = xp(1) * xp(1) + xp(3) * xp(3);
    
    rp = sqrt(rp2);
    vrp = (xp(1) * xp(2) + xp(3) * xp(4)) / rp;
    thetap = atan2(xp(3), xp(1));
    
    dr = r - rp;
    dvr = vr - vrp;
    dtheta = wrap_minus_plus_pi(theta - thetap);
    dz = [ dr ; dvr ; dtheta ];
    
    H = [ xp(1) / rp    , 0          , xp(3) / rp  , 0          ;
          0             , xp(1) / rp , 0           , xp(3) / rp ;
          - xp(3) / rp2 , 0          , xp(1) / rp2 , 0          ];
    
    % S = H * Pp * H' + R;
    S = quad_mult(H, Pp) + R;
    % K = Pp * H' * inv(S);
    K = Pp * H' * inv_sym_3x3(S);
    
    x = xp + K * dz;
    % P = Pp - K * S * K';
    P = Pp - quad_mult(K, S);
    
    s_k.t = t;
    s_k.x = x;
    s_k.P = P;
    
end
