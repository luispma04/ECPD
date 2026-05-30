function [u0, exitflag] = mpc_solve_sparse(x0, H, R, A, B, C, lb, ub, y_max_inc)
%% Sparse MPC Solver (Quadprog Formulation)
%--------------------------------------------------------------------------
% Solves the finite-horizon quadratic MPC regulator problem using a sparse
% formulation. Keeps both state and control variables as explicit 
% optimization variables with dynamics enforced as equality constraints.
%
% Supports bounded control inputs and constrained outputs.
%
% Inputs:
%   x0          - Current state vector (n x 1)
%   H           - Prediction horizon (positive integer)
%   R           - Control weight (positive scalar)
%   A           - State-transition matrix (n x n)
%   B           - Input matrix (n x 1)
%   C           - Output matrix (1 x n)
%   lb          - Lower bound on control sequence U (H x 1) (Optional)
%   ub          - Upper bound on control sequence U (H x 1) (Optional)
%   y_max_inc   - Maximum allowed incremental output (scalar) (Optional)
%
% Outputs:
%   u0          - First optimal control action to apply (scalar)
%   exitflag    - Solver exit flag from quadprog
%
% Matlab Toolboxes: Optimization Toolbox
% Functions: quadprog
%--------------------------------------------------------------------------

%% 1. Default Arguments

% If optional arguments are not provided, define them as empty
if nargin < 7
    lb = []; 
end
if nargin < 8
    ub = []; 
end
if nargin < 9
    y_max_inc = []; % Empty means no output constraint
end

use_output_constraint = ~isempty(y_max_inc);

%% 2. Dimensionality Setup

n   = size(A, 1);    % State vector dimension
N_x = (H + 1) * n;   % Number of state optimization variables
N_u = H;             % Number of control optimization variables
N_z = N_x + N_u;     % Total decision variables in the quadprog vector

%% 3. Cost Formulation

% Stage output cost: Q_stage = C'*C (n x n positive semidefinite)
Q_stage = C' * C;

% Qtilde = blkdiag(0_n, Q_stage, ..., Q_stage) with (H+1) diagonal blocks.
% The first block is 0 because xhat(0) does not appear in the cost.
Qtilde = blkdiag(zeros(n), kron(eye(H), Q_stage));

% Control penalty matrix
Rtilde = R * eye(H);

% Quadprog minimizes (1/2)*z'*F*z + f'*z
% From J = X'*Qtilde*X + U'*Rtilde*U
F = 2 * blkdiag(Qtilde, Rtilde);
f = zeros(N_z, 1);

%% 4. Equality Constraints (Dynamics & Initial Condition)
% z is structured as: [xhat(0); xhat(1); ...; xhat(H); uhat(0); ...; uhat(H-1)]

N_eq = (H + 1) * n; % Total number of equality constraint rows
Aeq  = zeros(N_eq, N_z);
beq  = zeros(N_eq, 1);

% Initial condition block: xhat(0) = x0
Aeq(1:n, 1:n) = eye(n);
beq(1:n)      = x0;

% Dynamics blocks: xhat(ii+1) - A*xhat(ii) - B*uhat(ii) = 0
for ii = 0:H-1
    r1 = n + ii*n + 1;       % First row for xhat(ii+1) constraint
    r2 = n + (ii+1)*n;       % Last row for xhat(ii+1) constraint

    % Coefficient of xhat(ii+1): +I
    c1_xi1 = (ii+1)*n + 1;
    c2_xi1 = (ii+2)*n;
    Aeq(r1:r2, c1_xi1:c2_xi1) = eye(n);

    % Coefficient of xhat(ii): -A
    c1_xi = ii*n + 1;
    c2_xi = (ii+1)*n;
    Aeq(r1:r2, c1_xi:c2_xi) = -A;

    % Coefficient of uhat(ii): -B (scalar control means a single column)
    c_ui = N_x + ii + 1;
    Aeq(r1:r2, c_ui) = -B;
end

%% 5. Control Variable Bounds

% State variables remain unconstrained; bounds apply only to U block
if ~isempty(lb) || ~isempty(ub)
    if isempty(lb)
        lb = -inf(N_u, 1); 
    end
    if isempty(ub)
        ub =  inf(N_u, 1); 
    end
    
    lb_z = [-inf(N_x, 1); lb];
    ub_z = [ inf(N_x, 1); ub];
else
    lb_z = [];
    ub_z = [];
end

%% 6. Output Constraints (Inequalities)

if use_output_constraint
    % C_til maps X to outputs [yhat(1); ...; yhat(H)]
    % xhat(0) does not appear in the cost -> prepend zero block
    C_til = kron(eye(H), C);          % H x H*n
    C_til = [zeros(H, n), C_til];     % H x (H+1)*n 

    % Full inequality: [C_til, 0_{H x H}] * z <= y_max_inc * ones(H,1)
    A_ineq = [C_til, zeros(H, N_u)];  % H x N_z
    b_ineq = y_max_inc * ones(H, 1);  % H x 1
else
    A_ineq = [];
    b_ineq = [];
end

%% 7. Solve Optimization Problem

opts = optimoptions('quadprog', 'Display', 'off');
[U_opt, ~, exitflag] = quadprog(F, f, A_ineq, b_ineq, Aeq, beq, lb_z, ub_z, [], opts);

%% 8. Receding Horizon Step

% U vector starts at index N_x + 1 in the optimized z vector
u0 = U_opt(N_x + 1);

end