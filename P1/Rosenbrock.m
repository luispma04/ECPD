function f = Rosenbrock(x)
% ROSENBROCK Computes the value of the Rosenbrock function in 2D.
%
% Input:
%   x - A 2-dimensional vector [x1; x2]
%
% Output:
%   f - Value of the function at point x
%
% Formula:
%   f(x1, x2) = 100*(x2 - x1^2)^2 + (1 - x1)^2
%
% The global unconstrained minimum is located at [1; 1] where f = 0.
%--------------------------------------------------------------------------

% Extract vector elements for clarity
x1 = x(1);
x2 = x(2);

% Compute Rosenbrock function value
f = 100 * (x2 - x1^2)^2 + (1 - x1)^2;

end