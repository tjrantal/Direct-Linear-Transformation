%solve equations
%-2x^2 + 3xy   + 4 sin(y) = 6
% 3x^2 - 2xy^2 + 3 cos(x) = -4
%use fsolve to solve the roots at 1,2
[x, fval, info] = fsolve (@fTest, [1; 2])

