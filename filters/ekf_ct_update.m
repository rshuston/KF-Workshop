# This EKF implements the "constant turn" target model with acceleration
# disturbances.
# This version of the EKF manually accounts for angles crossing the
# [-pi, pi) boundary.

function s_k = ekf_ct_update(t, s_km1, r, theta, proc_vars, meas_vars)
    
    var_a = proc_vars(1);       # cartesian acceleration disturbance
    var_alpha = proc_vars(2);   # angular acceleration disturbance
    var_r = meas_vars(1);
    var_theta = meas_vars(2);
    
    T = t - s_km1.t;
    T2 = T * T;
    
    x = s_km1.x;
    P = s_km1.P;
    
    Q = [var_a , 0 ,     0;
         0 ,     var_a , 0;
         0 ,     0 ,     var_alpha];
    
    R = [ var_r , 0         ;
          0     , var_theta ];
    
    # Prediction
    
    w = x(5);
    if (w == 0)
        cos_wT = 1;
        sin_wT = 0;
        sin_wT_div_w = T;
        one_m_cos_wT_div_w = 0;
        dndw = -0.5 * T2 * x(4);
        dvndw = -T * x(4);
        dedw = 0.5 * T2 * x(2);
        dvedw = T * x(2);
    else
        cos_wT = cos(w * T);
        sin_wT = sin(w * T);
        sin_wT_div_w = sin_wT / w;
        one_m_cos_wT_div_w = (1 - cos_wT) / w;
        temp1 = (w * T * cos_wT - sin_wT) / (w * w);
        temp2 = (w * T * sin_wT -1.0 + cos_wT) / (w * w);
        vn = x(2);
        ve = x(4);
        dndw = vn * temp1 - ve * temp2;
        dvndw = -vn * T * sin_wT - ve * T * cos_wT;
        dedw = vn * temp2 + ve * temp1;
        dvedw = vn * T * cos_wT - ve * T * sin_wT;
    endif
    
    Phi = [1 , sin_wT_div_w ,       0 , -one_m_cos_wT_div_w , 0;
           0 , cos_wT ,             0 , -sin_wT ,             0;
           0 , one_m_cos_wT_div_w , 1 , sin_wT_div_w ,        0;
           0 , sin_wT ,             0 , cos_wT ,              0;
           0 , 0 ,                  0 , 0 ,                   1];
    
    F = Phi;
    F(1,5) = dndw;
    F(2,5) = dvndw;
    F(3,5) = dedw;
    F(4,5) = dvedw;
    
    Gamma = [0.5 * T2 , 0 ,        0;
             T ,        0 ,        0;
             0 ,        0.5 * T2 , 0;
             0 ,        T ,        0;
             0 ,        0 ,        T];
    
    xp = Phi * x;
    Pp = F * P * F' + Gamma * Q * Gamma';
    
    n2e2 = xp(1) * xp(1) + xp(3) * xp(3);
    
    rp = sqrt(n2e2);
    thetap = atan2(xp(3), xp(1));
    
    dr = r - rp;
    dtheta = wrap_minus_plus_pi(theta - thetap);
    dz = [ dr ; dtheta ];
    
    H = [ xp(1) / sqrt(n2e2) , 0 , xp(3) / sqrt(n2e2) , 0 , 0 ;
          - xp(3) / n2e2 ,     0 , xp(1) / n2e2       , 0 , 0 ];
    
    S = H * Pp * H' + R;
    K = Pp * H' * inv(S);
    
    x = xp + K * dz;
    P = Pp - K * S * K';
    
    s_k.t = t;
    s_k.x = x;
    s_k.P = P;
    
endfunction
