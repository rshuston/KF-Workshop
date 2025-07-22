% This version of the UKF implements the classic Julier and Uhlmann algorithm.
% It uses direction cosine angle measurements so that we can easily handle
% angles crossing the [-pi, pi) boundary.

function s_k = ukf_ju_cvdc_update(t, s_km1, r, theta, proc_vars, meas_vars)
    
    var_p = proc_vars(1);
    var_r = meas_vars(1);
    var_dc = meas_vars(2);
    
    T = t - s_km1.t;
    
    x = s_km1.x;
    P = s_km1.P;
    
    N = s_km1.N;
    M = s_km1.M;
    N_plus_kappa = s_km1.N_plus_kappa;
    W1 = s_km1.W1;
    Wj = s_km1.Wj;
    
    two_N_plus_1 = 2 * N + 1;
    
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
    
    % Compute sigma points
            
    L = chol(N_plus_kappa * P, "lower");
    
    X = zeros(N, two_N_plus_1);
    X(:,1) = x;
    for i = 2:N+1
        X(:,i)     = x + L(:,i-1);
        X(:,i + N) = x - L(:,i-1);
    end
    
    % Compute predicted state and state covariance, saving intermediate results
    
    f_X = zeros(N,two_N_plus_1);
    for j = 1:two_N_plus_1
        f_X(:,j) = Phi * X(:,j);
    end
    
    xp = zeros(N,1);
    for j = 1:two_N_plus_1
        if (j == 1)
            Wm = W1;
        else
            Wm = Wj;
        end
        xp = xp + Wm * f_X(:,j);
    end
    
    Pp = zeros(N,N);
    for j = 1:two_N_plus_1
        if (j == 1)
            Wc = W1;
        else
            Wc = Wj;
        end
        dx = f_X(:,j) - xp;
        Pp = Pp + Wc * (dx * dx');
    end
    Pp = Pp + Q;
    
    % Compute predicted measurement, saving intermediate results
    
    zp = zeros(M,1);
    h_X = zeros(M, two_N_plus_1);
    for j = 1:two_N_plus_1
        f_X_j = f_X(:,j);
        if (j == 1)
            Wm = W1;
        else
            Wm = Wj;
        end
        rho_j = sqrt(f_X_j(1) * f_X_j(1) + f_X_j(3) * f_X_j(3));
        h_X_j = [ rho_j          ;
                  f_X_j(1) / rho_j ;
                  f_X_j(3) / rho_j ];
        h_X(:,j) = h_X_j;
        zp = zp + Wm * h_X_j;
    end
    
    % Normalize mean direction cosines adjusted from weighting
    c = zp(2);
    s = zp(3);
    d = sqrt(c * c + s * s);
    zp(2) = c / d;
    zp(3) = s / d;
    
    % Compute covariance matrices
    
    Pzz = zeros(M,M);
    Pxz = zeros(N,M);
    for j = 1:two_N_plus_1
        f_X_j = f_X(:,j);            
        h_X_j = h_X(:,j);
        if (j == 1)
            Wc = W1;
        else
            Wc = Wj;
        end

        dz = h_X_j - zp;
        dx = f_X_j - xp;
        
        Pzz = Pzz + Wc * (dz * dz');
        Pxz = Pxz + Wc * (dx * dz');
    end
    
    % Correction
    
    S = Pzz + R;
    K = Pxz * inv(S);
    
    dz = z - zp;
    
    x = xp + K * dz;
    P = Pp - K * S * K';

    % Simple hack to maintain symmetry
    P = 0.5 * (P + P');
    
    s_k.t = t;
    s_k.x = x;
    s_k.P = P;
    
    % Carry these forward to the next state update
    s_k.N = N;
    s_k.M = M;
    s_k.N_plus_kappa = N_plus_kappa;
    s_k.W1 = W1;
    s_k.Wj = Wj;
    
end
