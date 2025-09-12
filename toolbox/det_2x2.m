%
% Determinant of 2x2 matrix ... direct computation
%

function det_A = det_2x2(A)
    
    % A =
    %    a11 a12
    %    a21 a22
    
    det_A = A(1,1) * A(2,2) - A(1,2) * A(2,1);
    
end
