%
% Inversion of symmetric 3x3 matrix
%

function inv_A = inv_sym_3x3(A)
    
    % A =
    %    a11 a12 a13
    %    a12 a22 a23
    %    a13 a23 a33
    
    a11 = A(1,1);
    a12 = A(1,2);
    a13 = A(1,3);
    
    a22 = A(2,2);
    a23 = A(2,3);
    
    a33 = A(3,3);
    
    % Minors
    
    m11 = a22 * a33 - a23 * a23;
    m12 = a12 * a33 - a13 * a23;
    m13 = a12 * a23 - a13 * a22;
    
    m22 = a11 * a33 - a13 * a13;
    m23 = a11 * a23 - a13 * a12;
    
    m33 = a11 * a22 - a12 * a12;
    
    % Cofactors
    
    c11 = m11;
    c12 = -m12;
    c13 = m13;
    
    c22 = m22;
    c23 = -m23;
    
    c33 = m33;
    
    % Determinant
    
    D = a11 * c11 + a12 * c12 + a13 * c13;
    
    % Adjoint scaled by 1/D
    
    inv_A = (1 / D) * [ c11 c12 c13 ;
                        c12 c22 c23 ;
                        c13 c23 c33 ];
    
end
