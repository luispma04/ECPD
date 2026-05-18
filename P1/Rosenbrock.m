function f = Rosenbrock(x)
%   This function computes the Rosenbrock function in 2 dimensions at a point.
%   The argument x of the function is a 2-dimensional vector;
%   The function computes the value of the Rosenbrock function at x.
%
%   f(x1, x2) = 100*(x2 - x1^2)^2 + (1 - x1)^2
%
%   The global (unconstrained) minimum is at x* = [1; 1] with f(x*) = 0.
%
%   IST - MEEC - Distributed Predictive Control and Estimation
%   Afonso Botelho, Joao Miranda Lemos, 2025
%--------------------------------------------------------------------------

f = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;

end
%--------------------------------------------------------------------------
% End of file
