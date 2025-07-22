%
% Quadratic Multiply - determines C = A B A', where B is a symmetric matrix,
% and the symmetry of C is preserved
%
% A: m x n
% B: n x n
% C: m x m
%

function C = quad_mult(A, B)
    sA = size(A);
    sB = size(B);
    C = 0;
    if (sA(2) == sB(1) && sB(1) == sB(2))
        m = sA(1);
        n = sA(2);
        C = zeros(m);
        for i = 1:m
            sigma = 0;
            for j = 1:n
                for k = 1:n
                    sigma = sigma + A(i,j) * B(j,k) * A(i,k);
                end
            end
            C(i,i) = sigma;
            for l = (i+1):m
                sigma = 0;
                for j = 1:n
                    for k = 1:n
                        sigma = sigma + A(i,j) * B(j,k) * A(l,k);
                    end
                end
                C(i,l) = sigma;
                C(l,i) = sigma;
            end
        end
    end
end
