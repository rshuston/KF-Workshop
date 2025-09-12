%
% Determinant of 3x3 matrix ... direct computation
%

function det_A = det_3x3(A)
    
    % A =
    %    a11 a12 a13
    %    a21 a22 a23
    %    a31 a32 a33
    
    m11 = A(2,2) * A(3,3) - A(2,3) * A(3,2);
    m12 = A(2,1) * A(3,3) - A(2,3) * A(3,1);
    m13 = A(2,1) * A(3,2) - A(2,2) * A(3,1);
    
    det_A = A(1,1) * m11 - A(1,2) * m12 + A(1,3) * m13;
    
end
