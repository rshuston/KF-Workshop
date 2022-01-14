function [n, vn, e, ve, Pne] = ekf_cv(t, r, theta, var_p, var_r, var_theta)
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
    
    for k = 2:length(t)
        printf("########## k = \%d\n", k);

        T = t(k) - t(k-1);
        
        Phi = [1 , T , 0 , 0;
               0 , 1 , 0 , 0;
               0 , 0 , 1 , T;
               0 , 0 , 0 , 1];
        
        Phi

        T2 = T * T;        
        T3d2 = T * T * T / 2;
        T4d4 = T * T * T * T / 4;
        
#        Q = [ T2 , T , T2 , T ;
#              T  , 1 , T  , 1 ;
#              T2 , T , T2 , T ;
#              T  , 1 , T  , 1 ] * var_p;
  
#        Q = [ T4d4 , T3d2 , T4d4 , T3d2 ;
#              T3d2 , T2   , T3d2 , T2   ;
#              T4d4 , T3d2 , T4d4 , T3d2 ;
#              T3d2 , T2   , T3d2 , T2   ] * var_p;
        
#        Q = [ T2 , T , 0  , 0 ;
#              T  , 1 , 0  , 0 ;
#              0  , 0 , T2 , T ;
#              0  , 0 , T  , 1 ] * var_p;
  
        Q = [ T4d4 , T3d2 , 0    , 0    ;
              T3d2 , T2   , 0    , 0    ;
              0    , 0    , T4d4 , T3d2 ;
              0    , 0    , T3d2 , T2   ] * var_p;
        
        Q
        
        R = [ var_r , 0         ;
              0     , var_theta ];
        
        R
               
        xp = Phi * x;
        Pp = Phi * P * Phi' + Q;
        
        xp
        Pp

        n2e2 = xp(1) * xp(1) + xp(3) * xp(3);
        
        rp = sqrt(n2e2);
        thetap = atan2(xp(3), xp(1));
        
        dr = r(k) - rp;
        dtheta = theta(k) - thetap;
        if (dtheta < -pi)
            dtheta += 2 * pi;
        elseif (dtheta > pi)
            dtheta -= 2 * pi;
        endif
        
        dz = [ dr ; dtheta ];
        
        dz
        
        H = [ xp(1) / sqrt(n2e2) , 0 , xp(3) / sqrt(n2e2) , 0 ;
              - xp(3) / n2e2 ,     0 , xp(1) / n2e2       , 0 ];
        
        H
        
        S = H * Pp * H' + R;
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
