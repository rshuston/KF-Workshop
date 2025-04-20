# This version of the UKF only approximates the measurement relation, and uses
# different weights for mean (Wm) and covariance (Wc) sums. It uses direction
# cosine angle measurements so that we can easily handle angles crossing the
# [-pi, pi) boundary.

function s_k = ukf_cvdc_init(t, r, theta, vars)
    
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
    two_N_plus_1 = 2 * N + 1;
    alpha = 0.5;
    beta = 2;
    kappa = 0;
    lambda = alpha * alpha * (N + kappa) - N;
    Wm0 = lambda / (N + lambda);
    Wc0 = Wm0 + 1 - alpha * alpha + beta;
    Wi = 1 / (2 * (N + lambda));
    sqrt_N_plus_lambda = sqrt(alpha * alpha * (N + kappa));
    
    Wm = zeros(1,two_N_plus_1);
    Wc = zeros(1,two_N_plus_1);
    Wm(1) = Wm0;
    Wc(1) = Wc0;
    for i = 2:two_N_plus_1
        Wm(i) = Wi;
        Wc(i) = Wi;
    endfor
    
    s_k.N = N;
    s_k.M = M;
    s_k.sqrt_N_plus_lambda = sqrt_N_plus_lambda;
    s_k.Wm = Wm;
    s_k.Wc = Wc;
    
endfunction
