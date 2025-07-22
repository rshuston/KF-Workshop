function inv_M = inv_sym_2x2(M)
    
    % M =
    %    m11 m12
    %    m12 m22
    
    m11 = M(1,1);
    m12 = M(1,2);
    
    m22 = M(2,2);
    
    D = m11 * m22 - m12 * m12;
    
    inv_M = (1 / D) * [  m22 -m12 ;
                        -m12  m11 ];
    
end
