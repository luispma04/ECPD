function [KRH, M, W, Pi] = compute_KRH(A, b, C, H, R)
% COMPUTE_KRH  Receding-Horizon state-feedback gain for a SISO plant.
%
%   [KRH, M, W, Pi] = compute_KRH(A, b, C, H, R)
%
%   Computes the optimal receding-horizon gain KRH for the scalar-input
%   plant
%       x(t+1) = A x(t) + b u(t),   y(t) = C x(t),
%   by minimising the finite-horizon quadratic cost
%       J = sum_{i=0}^{H-1} [ y(t+i+1)^2 + R*u(t+i)^2 ]
%   over the control sequence U = [u(t); ...; u(t+H-1)].
%
%   The optimal unconstrained solution is
%       U* = -M^{-1} W' Pi x(t),
%   and only the first element is applied (receding horizon):
%       u(t) = -KRH x(t),   KRH = e1 * M^{-1} * W' * Pi.
%
%   Inputs
%     A  – n×n state-transition matrix
%     b  – n×1 input matrix  (scalar input assumed)
%     C  – 1×n output matrix (scalar output assumed)
%     H  – prediction horizon (positive integer)
%     R  – positive scalar control weight
%
%   Outputs
%     KRH – 1×n receding-horizon feedback gain
%     M   – H×H Hessian matrix  M = W'W + R*I
%     W   – H×H lower-triangular prediction matrix
%     Pi  – H×n initial-state prediction matrix

n = size(A, 1);

%% Build W  (H×H lower-triangular Toeplitz in C*A^k*b)
%  W(i,j) = C * A^(i-j) * b   for  i >= j,  else 0
W  = zeros(H, H);
for ii = 1:H
    for jj = 1:ii
        W(ii, jj) = C * (A^(ii-jj)) * b;
    end
end

%% Build Pi  (H×n)
%  Pi(i,:) = C * A^i
Pi = zeros(H, n);
for ii = 1:H
    Pi(ii, :) = C * (A^ii);
end

%% Hessian
M = W' * W + R * eye(H);

%% Receding-horizon gain  (only the first row of M^{-1}W'Pi)
e1  = zeros(1, H);
e1(1) = 1;                    % selects u(t) from U*
KRH = e1 * (M \ (W' * Pi));  % 1×n row vector

end
