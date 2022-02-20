#
# Wraps X into the [-pi,  pi) interval
#

function xwrap = wrapToPi(x)
    xwrap = mod(x,2*pi);
    idx = xwrap >= pi;
    xwrap(idx) -= 2*pi;
endfunction

#
# Example:
#
# wrapToPi ([-3*pi, -pi, -pi-1, 0; pi-1, pi, pi+1, 3*pi])
# ans =
#
#   -3.1416  -3.1416   2.1416        0
#    2.1416   3.1416  -2.1416   3.1416
#
