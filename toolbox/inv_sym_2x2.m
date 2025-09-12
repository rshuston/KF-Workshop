%
% Inversion of symmetric 2x2 matrix
%

function inv_A = inv_sym_2x2(A)
    
    % A =
    %    a11 a12
    %    a12 a22
    
    a11 = A(1,1);
    a12 = A(1,2);
    
    a22 = A(2,2);
    
    D = a11 * a22 - a12 * a12;
    
    inv_A = (1 / D) * [  a22 -a12 ;
                        -a12  a11 ];
    
end
