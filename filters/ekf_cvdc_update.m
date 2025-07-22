% This version of the EKF uses direction cosine angle measurements so that we
% can easily handle angles crossing the [-pi, pi) boundary.

function s_k = ekf_cvdc_update(t, s_km1, r, theta, proc_vars, meas_vars)
    
    var_p = proc_vars(1);
    var_r = meas_vars(1);
    var_dc = meas_vars(2);
    
    T = t - s_km1.t;
    
    x = s_km1.x;
    P = s_km1.P;
    
    z = [ r      ;
          cos(theta) ;
          sin(theta) ];
    
    % Heuristic covariance, but numerically behaved
    R = [ var_r , 0      , 0      ;
          0     , var_dc , 0      ;
          0     , 0      , var_dc ];
    
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
    
    xp = Phi * x;
    % Pp = Phi * P * Phi' + Q;
    Pp = quad_mult(Phi, P) + Q;
    
    rp = sqrt(xp(1) * xp(1) + xp(3) * xp(3));
    
    zp = [ rp         ;
           xp(1) / rp ;
           xp(3) / rp ];
           
    % Normalize direction cosines adjusted from roundoff
    c = zp(2);
    s = zp(3);
    d = sqrt(c * c + s * s);
    zp(2) = c / d;
    zp(3) = s / d;
    
    dz = z - zp;
    
    rp3 = rp * rp * rp;
    H = [ xp(1) / rp            , 0 , xp(3) / rp            , 0 ;
          xp(3) * xp(3) / rp3   , 0 , - xp(1) * xp(3) / rp3 , 0 ;
          - xp(1) * xp(3) / rp3 , 0 , xp(1) * xp(1) / rp3   , 0 ];
    
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
