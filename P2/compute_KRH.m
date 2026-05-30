function [KRH, M, W, Pi] = compute_KRH(A, b, C, H, R)
% COMPUTE_KRH  Receding-Horizon state-feedback gain for a SISO plant.
%
%   Computes the optimal receding-horizon gain KRH for the system:
%       x(t+1) = A x(t) + b u(t)
%       y(t)   = C x(t)
%
%   Minimizes the finite-horizon quadratic cost:
%       J = sum_{ii=0}^{H-1} [ y(t+ii+1)^2 + R*u(t+ii)^2 ]
%
%   Inputs:
%     A  - State-transition matrix (n x n)
%     b  - Input matrix (n x 1)
%     C  - Output matrix (1 x n)
%     H  - Prediction horizon (positive integer)
%     R  - Control weight (positive scalar)
%
%   Outputs:
%     KRH - Receding-horizon feedback gain (1 x n)
%     M   - Hessian matrix (H x H)
%     W   - Lower-triangular prediction matrix (H x H)
%     Pi  - Initial-state prediction matrix (H x n)
%--------------------------------------------------------------------------

% Get state dimension
n = size(A, 1);

%% Build W: Lower-triangular Toeplitz matrix
% W(ii,jj) = C * A^(ii-jj) * b for ii >= jj, else 0
W = zeros(H, H);
for ii = 1:H
    for jj = 1:ii
        W(ii, jj) = C * (A^(ii-jj)) * b;
    end
end

%% Build Pi: Initial-state prediction matrix
% Pi(ii,:) = C * A^ii
Pi = zeros(H, n);
for ii = 1:H
    Pi(ii, :) = C * (A^ii);
end

%% Hessian Matrix Computation
% M = W'W + R*I
M = W' * W + R * eye(H);

%% Receding-Horizon Gain Calculation
% We apply only the first element of the control sequence U*
% U* = -M^{-1} * W' * Pi * x(t)
% KRH = e1 * M^{-1} * W' * Pi, where e1 = [1 0 ... 0]
e1 = zeros(1, H);
e1(1) = 1;

% KRH is computed using backslash operator for numerical stability
KRH = e1 * (M \ (W' * Pi));

end