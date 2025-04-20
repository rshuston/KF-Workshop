# This version of the UKF implements the classic Julier and Uhlmann algorithm.
# It uses direction cosine angle measurements so that we can easily handle
# angles crossing the [-pi, pi) boundary.

function s_k = ukf_ju_cvdc_init(t, r, theta, vars)
    
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
    
    N = 4; # states
    M = 3; # measurements
    N_plus_kappa = 3;
    kappa = N_plus_kappa - N;
    W1 = kappa / N_plus_kappa;
    Wj = 1 / (2 * N_plus_kappa);
    
    s_k.N = N;
    s_k.M = M;
    s_k.N_plus_kappa = N_plus_kappa;
    s_k.W1 = W1;
    s_k.Wj = Wj;
    
endfunction
