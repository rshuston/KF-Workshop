# Standard linear Kalman filter ... nothing fancy

function s_k = kf_cv_init(t, r, theta, vars)
    
    n = r * cos(theta);
    vn = 0;
    e = r * sin(theta);
    ve = 0;
    
    var_ne = vars(1);
    var_vnve = vars(2);
    P0 = [ var_ne , 0 ,        0 ,      0 ;
           0 ,      var_vnve , 0 ,      0 ;
           0 ,      0 ,        var_ne , 0 ;
           0 ,      0 ,        0 ,      var_vnve ];
    
    s_k.t = t;
    s_k.x = [n; vn; e; ve];
    s_k.P = P0;
    
endfunction
