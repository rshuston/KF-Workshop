# This EKF implements the "constant turn" target model with acceleration
# disturbances.
# This version of the EKF uses direction cosine angle measurements so that we
# can easily handle angles crossing the [-pi, pi) boundary.

function s_k = ekf_ctdc_init(t, r, theta, vars)
    
    n = r * cos(theta);
    vn = 0;
    e = r * sin(theta);
    ve = 0;
    w = 0;
    
    var_ne = vars(1);
    var_vnve = vars(2);
    var_w = vars(3);
    P0 = [ var_ne , 0 ,        0 ,      0 ,        0 ;
           0 ,      var_vnve , 0 ,      0 ,        0 ;
           0 ,      0 ,        var_ne , 0 ,        0 ;
           0 ,      0 ,        0 ,      var_vnve , 0 ;
           0 ,      0 ,        0 ,      0 ,        var_w ];
    
    s_k.t = t;
    s_k.x = [n; vn; e; ve; w];
    s_k.P = P0;
    
endfunction
