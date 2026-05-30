function u0 = mpc_solve(x0, H, R, A, B, C, lb, ub)
%% Dense MPC Solver (Quadprog Formulation)
%--------------------------------------------------------------------------
% Solves the finite-horizon quadratic MPC regulator problem (dense 
% formulation) using quadratic programming (quadprog).
%
% Minimizes the cost function:
%   J = Y'Y + R*U'U
%
% Where Y = W*U + Pi*x0. This expands to the quadprog standard form:
%   min_U (1/2)*U'*F*U + f'*U
%
% Inputs:
%   x0 - Current state vector (n x 1)
%   H  - Prediction horizon (positive integer)
%   R  - Control weight (positive scalar)
%   A  - State-transition matrix (n x n)
%   B  - Input matrix (n x 1)
%   C  - Output matrix (1 x n)
%   lb - Lower bound on control sequence U (H x 1) (Optional)
%   ub - Upper bound on control sequence U (H x 1) (Optional)
%
% Outputs:
%   u0 - First optimal control action to apply (scalar)
%
% Note: This function supports both unconstrained (by omitting lb/ub) and 
% constrained formulations.
%
% Matlab Toolboxes: Optimization Toolbox
% Functions: quadprog
%--------------------------------------------------------------------------

%% 1. Default Arguments (Unconstrained default)

% If boundaries are not provided, define them as empty arrays
if nargin < 7
    lb = []; 
end
if nargin < 8
    ub = []; 
end

%% 2. Build Prediction Matrices

n = size(A, 1); % State vector dimension

% Initialize prediction matrices
W  = zeros(H, H);
Pi = zeros(H, n);

% Build W: Lower-triangular matrix where W(ii,jj) = C * A^(ii-jj) * B
for ii = 1:H
    for jj = 1:ii
        W(ii, jj) = C * (A^(ii-jj)) * B;
    end
end

% Build Pi: Initial state prediction matrix where Pi(ii,:) = C * A^ii
for ii = 1:H
    Pi(ii, :) = C * (A^ii);
end

%% 3. Formulate Quadprog Cost Matrices

% Hessian of the cost w.r.t U: M = W'*W + R*I
M = W' * W + R * eye(H);

% Quadprog minimizes: (1/2)*x'*H_qp*x + f_qp'*x
% From J = U'*M*U + 2*x0'*Pi'*W*U + const
% We match forms by setting F = 2*M and f = 2*W'*Pi*x0
F = 2 * M;
f = 2 * W' * Pi * x0;

%% 4. Solve Optimization Problem

% Suppress quadprog console output for cleaner loop execution
opts = optimoptions('quadprog', 'Display', 'off');

% Solve for the optimal control sequence U*
U_opt = quadprog(F, f, [], [], [], [], lb, ub, [], opts);

%% 5. Receding Horizon Step

% Apply only the first computed control action
u0 = U_opt(1);

end