% Standard linear Kalman filter ... nothing fancy

function s_k = kf_cv_update(t, s_km1, r, theta, proc_vars, meas_vars)
    
    var_p = proc_vars(1);
    var_r = meas_vars(1);
    var_theta = meas_vars(2);
    
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
    
    xp = Phi * x;
    % Pp = Phi * P * Phi' + Q;
    Pp = quad_mult(Phi, P) + Q;
    
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    
    zn = r * cos_theta;
    ze = r * sin_theta;
    
    dz = [ zn - xp(1) ;
           ze - xp(3) ];
    
    M = [ cos_theta, - r * sin_theta ;
          sin_theta, r * cos_theta   ];
    
    R = [ var_r , 0         ;
          0     , var_theta ];
    
    % Rz = M * R * M';
    Rz = quad_mult(M, R);
    
    H = [ 1 , 0 , 0 , 0 ;
          0 , 0 , 1 , 0 ];
    
    % S = H * Pp * H' + Rz;
    S = quad_mult(H, Pp) + Rz;
    % K = Pp * H' * inv(S);
    K = Pp * H' * inv_sym_2x2(S);
    
    x = xp + K * dz;
    % P = Pp - K * S * K';
    P = Pp - quad_mult(K, S);
    
    s_k.t = t;
    s_k.x = x;
    s_k.P = P;
    
end
