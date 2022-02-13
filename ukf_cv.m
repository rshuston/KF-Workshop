function [n, vn, e, ve, Pne] = ukf_cv(t, r, theta, var_p, var_r, var_theta)
    n  = zeros(1, length(t));
    vn = zeros(1, length(t));
    e  = zeros(1, length(t));
    ve = zeros(1, length(t));
    Pne = cell(1, length(t));

    printf("########## k = \%d\n", 1);
    
    x = [ r(1) * cos(theta(1)) ;
          0                    ;
          r(1) * sin(theta(1)) ;
          0                    ];
    P = [1 , 0 , 0 , 0;
         0 , 1 , 0 , 0;
         0 , 0 , 1 , 0;
         0 , 0 , 0 , 1];
    
    x
    P
    
    n(1)  = x(1);
    vn(1) = x(2);
    e(1)  = x(3);
    ve(1) = x(4);
    
    Pne{1} = P;
    
    N = 4;
    M = 2;
    two_N_plus_1 = 2 * N + 1;
    alpha = 0.5;
    beta = 2;
    kappa = 0;
    lambda = alpha * alpha * (N + kappa) - N;
    Wm0 = lambda / (N + lambda);
    Wc0 = Wm0 + 1 - alpha * alpha + beta;
    Wi = 1 / (2 * (N + lambda));
    sqrt_N_plus_lambda = sqrt(alpha * alpha * (N + kappa));
    
    N
    M
    two_N_plus_1
    lambda
    sqrt_N_plus_lambda
    Wm0
    Wc0
    
    R = [ var_r , 0         ;
          0     , var_theta ];
    
    for k = 2:length(t)
        printf("########## k = \%d\n", k);
        
        # Prediction

        T = t(k) - t(k-1);
        
        Phi = [1 , T , 0 , 0;
               0 , 1 , 0 , 0;
               0 , 0 , 1 , T;
               0 , 0 , 0 , 1];
        
        Phi

        T2 = T * T;        
        T3d2 = T * T * T / 2;
        T4d4 = T * T * T * T / 4;
        
#        Q = [ T2 , T , 0  , 0 ;
#              T  , 1 , 0  , 0 ;
#              0  , 0 , T2 , T ;
#              0  , 0 , T  , 1 ] * var_p;
  
        Q = [ T4d4 , T3d2 , 0    , 0    ;
              T3d2 , T2   , 0    , 0    ;
              0    , 0    , T4d4 , T3d2 ;
              0    , 0    , T3d2 , T2   ] * var_p;
        
        Q
               
        xp = Phi * x;
        Pp = Phi * P * Phi' + Q;
        
        xp
        Pp
        
        # Compute sigma points and weights
        
        L = sqrt_N_plus_lambda * chol(Pp, "lower");
        L
        
        X = zeros(N, two_N_plus_1);
        X(:,1) = xp;
        for i = 2:N+1
            X(:,i) = xp + L(:,i-1);
            X(:,i + N) = xp - L(:,i-1);
        endfor
        
        X
        
        Wm = zeros(1,two_N_plus_1);
        Wc = zeros(1,two_N_plus_1);
        Wm(1) = Wm0;
        Wc(1) = Wc0;
        for i = 2:two_N_plus_1
            Wm(i) = Wi;
            Wc(i) = Wi;
        endfor
        
        Wm
        Wc
        
        # Compute predicted measurement, saving intermediate results
        
        zp = zeros(M,1);
        h_X = zeros(M, two_N_plus_1);
        for i = 1:two_N_plus_1
            X_i = X(:,i);
            Wm_i = Wm(i);
            h_X_i = [ sqrt(X_i(1) * X_i(1) + X_i(3) * X_i(3)) ;
                      atan2(X_i(3), X_i(1))                   ];
            h_X(:,i) = h_X_i;
            zp += Wm_i * h_X_i;
        endfor
        
        # Normalize angle
        zp(2) = wrapToPi(zp(2));
        
        h_X
        zp
        
        # Compute covariance matrices
        
        Pxz = zeros(N,M);
        Pzz = zeros(M,M);
        for i = 1:two_N_plus_1
            X_i = X(:,i);
            h_X_i = h_X(:,i);
            Wc_i = Wc(i);
            
            dx = X_i - xp;
            dz = h_X_i - zp;
            dz(2) = wrapToPi(dz(2));
            
            Pxz += Wc_i * (dx * dz');
            Pzz += Wc_i * (dz * dz');
        endfor
        
        Pxz
        Pzz
        
        # Correction
        
        S = Pzz + R;
        K = Pxz * inv(S);
        
        S
        K
        
        dr = r(k) - zp(1);
        dtheta = wrapToPi(theta(k) - zp(2));
        dz = [ dr ; dtheta ];
        
        dz
        
        x = xp + K * dz;
        P = Pp - K * S * K';
        
        x
        P
        
        n(k)  = x(1);
        vn(k) = x(2);
        e(k)  = x(3);
        ve(k) = x(4);
        
        Pne{k} = P;
    endfor
endfunction
