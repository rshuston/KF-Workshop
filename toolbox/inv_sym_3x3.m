function inv_M = inv_sym_3x3(M)
    
    % M =
    %    m11 m12 m13
    %    m12 m22 m23
    %    m13 m23 m33
    
    m11 = M(1,1);
    m12 = M(1,2);
    m13 = M(1,3);
    
    m22 = M(2,2);
    m23 = M(2,3);
    
    m33 = M(3,3);
    
    % Cofactors
    
    c11 = m22 * m33 - m23 * m23;
    c12 = -(m12 * m33 - m13 * m23);
    c13 = m12 * m23 - m13 * m22;
    
    c22 = m11 * m33 - m13 * m13;
    c23 = -(m11 * m23 - m13 * m12);
    
    c33 = m11 * m22 - m12 * m12;
    
    % Determinant
    
    %D = m11 * (m22 * m33 - m23 * m23) - m12 * (m12 * m33 - m13 * m23) + m13 * (m12 * m23 - m13 * m22);
    D = m11 * c11 + m12 * c12 + m13 * c13;
    
    % Adjoint scaled by 1/D
    
    inv_M = (1 / D) * [ c11 c12 c13 ;
                        c12 c22 c23 ;
                        c13 c23 c33 ];
    
end
