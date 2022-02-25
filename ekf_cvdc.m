function [n, vn, e, ve, Pne] = ekf_cvdc(t, r, theta, P0, var_p, var_r, var_dc)
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
               
        xp = Phi * x;
        Pp = Phi * P * Phi' + Q;
        
        xp
        Pp

        rp = sqrt(xp(1) * xp(1) + xp(3) * xp(3));
        
        zp = [ rp         ;
               xp(1) / rp ;
               xp(3) / rp ];
               
        # Normalize direction cosines adjusted from roundoff
        c = zp(2);
        s = zp(3);
        d = sqrt(c * c + s * s);
        zp(2) = c / d;
        zp(3) = s / d;
        
        zp
        
        dz = z_k - zp;
        
        dz
        
        rp3 = rp * rp * rp;
        H = [ xp(1) / rp            , 0 , xp(3) / rp            , 0 ;
              xp(3) * xp(3) / rp3   , 0 , - xp(1) * xp(3) / rp3 , 0 ;
              - xp(1) * xp(3) / rp3 , 0 , xp(1) * xp(1) / rp3   , 0 ];
        
        H
        
        S = H * Pp * H' + R_k;
        K = Pp * H' * inv(S);
        
        S
        K
        
        x = xp + K * dz;
        P = (eye(4) - K * H) * Pp;
        
        x
        P
        
        n(k)  = x(1);
        vn(k) = x(2);
        e(k)  = x(3);
        ve(k) = x(4);
        
        Pne{k} = P;
    endfor
endfunction
