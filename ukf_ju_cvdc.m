# This version of the UKF implements the classic Julier and Uhlmann algorithm.
# It uses direction cosine angle measurements so that we can easily handle
# angles crossing the [-pi, pi) boundary.

function [n, vn, e, ve, Pne] = ukf_ju_cvdc(t, r, theta, P0, var_p, var_r, var_dc)
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
    P = P0;
    
    x
    P
    
    n(1)  = x(1);
    vn(1) = x(2);
    e(1)  = x(3);
    ve(1) = x(4);
    
    Pne{1} = P;
    
    N = 4;
    M = 3;
    two_N_plus_1 = 2 * N + 1;
        
    N
    M
    two_N_plus_1
        
    N_plus_kappa = 3;
    kappa = N_plus_kappa - N;
    W1 = kappa / N_plus_kappa;
    Wj = 1 / (2 * N_plus_kappa);

    N_plus_kappa
    kappa
    W1
    Wj
        
    for k = 2:length(t)
        printf("########## k = \%d\n", k);
        
        z_k = [ r(k)      ;
                cos(theta(k)) ;
                sin(theta(k)) ];
        
        # Heuristic covariance, but numerically behaved
        R_k = [ var_r , 0      , 0      ;
                0     , var_dc , 0      ;
                0     , 0      , var_dc ];
        
        z_k
        R_k
                
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
        
        # Compute sigma points
                
        L = chol(N_plus_kappa * P, "lower");
        
        L
        
        X = zeros(N, two_N_plus_1);
        X(:,1) = x;
        for i = 2:N+1
            X(:,i)     = x + L(:,i-1);
            X(:,i + N) = x - L(:,i-1);
        endfor
        
        X
        
        # Compute predicted state and state covariance, saving intermediate results
        
        f_X = zeros(N,two_N_plus_1);
        for j = 1:two_N_plus_1
            f_X(:,j) = Phi * X(:,j);
        endfor
        
        f_X
        
        xp = zeros(N,1);
        for j = 1:two_N_plus_1
            if (j == 1)
                Wm = W1
            else
                Wm = Wj;
            endif
            xp += Wm * f_X(:,j);
        endfor
        
        xp
        
        Pp = zeros(N,N);
        for j = 1:two_N_plus_1
            if (j == 1)
                Wc = W1
            else
                Wc = Wj;
            endif
            dx = f_X(:,j) - xp;
            Pp += Wc * (dx * dx');
        endfor
        Pp += Q;
        
        Pp
        
        # Compute predicted measurement, saving intermediate results
        
        zp = zeros(M,1);
        h_X = zeros(M, two_N_plus_1);
        for j = 1:two_N_plus_1
            f_X_j = f_X(:,j);
            if (j == 1)
                Wm = W1
            else
                Wm = Wj;
            endif
            rho_j = sqrt(f_X_j(1) * f_X_j(1) + f_X_j(3) * f_X_j(3));
            h_X_j = [ rho_j          ;
                      f_X_j(1) / rho_j ;
                      f_X_j(3) / rho_j ];
            h_X(:,j) = h_X_j;
            zp += Wm * h_X_j;
        endfor
        
        # Normalize mean direction cosines adjusted from weighting
        c = zp(2);
        s = zp(3);
        d = sqrt(c * c + s * s);
        zp(2) = c / d;
        zp(3) = s / d;
        
        h_X
        zp
        
        # Compute covariance matrices
        
        Pzz = zeros(M,M);
        Pxz = zeros(N,M);
        for j = 1:two_N_plus_1
            f_X_j = f_X(:,j);            
            h_X_j = h_X(:,j);
            if (j == 1)
                Wc = W1
            else
                Wc = Wj;
            endif

            dz = h_X_j - zp;
            dx = f_X_j - xp;
            
            Pzz += Wc * (dz * dz');
            Pxz += Wc * (dx * dz');
        endfor
        
        Pzz
        Pxz
        
        # Correction
        
        S = Pzz + R_k;
        K = Pxz * inv(S);
        
        S
        K
        
        dz = z_k - zp;
        
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
