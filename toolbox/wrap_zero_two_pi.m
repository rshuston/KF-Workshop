%
% Wraps X into the [0, 2 pi) interval
%

function xwrap = wrap_zero_two_pi(x)
    xwrap = mod(x, 2*pi);
end

%
% Example:
%
% wrap_zero_two_pi([-3*pi, -pi, -pi-1, 0; pi-1, pi, pi+1, 3*pi])
% ans =
%
%    3.1416   3.1416   2.1416        0
%    2.1416   3.1416   2.1416   3.1416
%
