function [u0, exitflag] = mpc_solve_sparse(x0, H, R, A, B, C, lb, ub, y_max_inc, const_type)
%% Regularized Sparse MPC Solver (Quadprog Formulation)
%--------------------------------------------------------------------------
% Solves the finite-horizon quadratic MPC regulator problem using a sparse
% formulation. Keeps state variables as explicit optimization variables 
% with dynamics enforced as equality constraints. Includes Tikhonov 
% regularization for numerical stability.
%
% Supports bounded control inputs and constrained outputs (both hard and 
% soft constraint formulations).
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
%   const_type  - Constraint type: 0 for Hard, 1 for Soft (Optional, def: 0)
%
% Outputs:
%   u0          - First optimal control action to apply (scalar)
%   exitflag    - Solver exit flag from quadprog
%
% Matlab Toolboxes: Optimization Toolbox
% Functions: quadprog
%--------------------------------------------------------------------------

%% 1. Default Arguments

% If optional arguments are not provided, define them as empty/default
if nargin < 7
    lb = []; 
end
if nargin < 8
    ub = []; 
end
if nargin < 9
    y_max_inc = []; 
end
if nargin < 10
    const_type = 0; % Default to hard constraints
end

use_output_constraint = ~isempty(y_max_inc);

%% 2. Dimensionality Setup

n        = size(A, 1);  % State vector dimension
N_x      = (H + 1) * n; % Number of state variables across horizon
N_u      = H;           % Number of control variables across horizon
N_eta    = H;           % Number of slack variables for soft constraints
N_z      = N_x + N_u;   % Total decision variables (standard)
N_z_soft = N_z + N_eta; % Total decision variables (with slack)

%% 3. Cost Formulation

% Base penalty matrices
Q_stage = C' * C;
Qtilde  = blkdiag(zeros(n), kron(eye(H), Q_stage));
Rtilde  = R * eye(H);

% Regularization parameters
epsilon = 1e-6; % Tikhonov regularization for state block
alpha   = 1e6;  % Large penalty weight for slack variables (soft constraints)

% F = 2 * blkdiag(Qtilde, Rtilde) + epsilon * I
% Q_stage has rank 1, so Qtilde is rank-deficient (many zero eigenvalues). 
% The epsilon*I regularization makes F strictly positive definite so the 
% KKT system is well-conditioned. We only regularize the state block.
F = 2 * blkdiag(Qtilde + epsilon * eye(N_x), Rtilde);
f = zeros(N_z, 1);

% Augment cost matrices for soft constraints
if const_type == 1
    % Extend F to include quadratic penalty on slack variables (eta)
    F = blkdiag(F, alpha * eye(N_eta));
    f = [f; zeros(N_eta, 1)];
end

%% 4. Equality Constraints (Dynamics & Initial Condition)

N_eq = N_x; % (H+1)*n equality constraints
Aeq  = zeros(N_eq, N_z);
beq  = zeros(N_eq, 1);

% Initial condition: xhat(0) = x0
Aeq(1:n, 1:n) = eye(n);
beq(1:n)      = x0;

% System Dynamics: xhat(ii+1) - A*xhat(ii) - B*uhat(ii) = 0
for ii = 0:H-1
    r1 = n + ii*n + 1;
    r2 = n + (ii+1)*n;
    
    Aeq(r1:r2, (ii+1)*n+1:(ii+2)*n) =  eye(n); % +I on xhat(ii+1)
    Aeq(r1:r2, ii*n+1:(ii+1)*n)     = -A;      % -A on xhat(ii)
    Aeq(r1:r2, N_x+ii+1)            = -B;      % -B on uhat(ii)
end

% Augment equality constraints for soft constraints
if const_type == 1
    % Extend Aeq with zero columns for the slack variables (eta)
    % The RHS vector (beq) remains unchanged
    Aeq = [Aeq, zeros(N_eq, N_eta)];
end

%% 5. Control and Slack Variable Bounds

if const_type == 0
    % Hard constraints boundary generation
    if ~isempty(lb) || ~isempty(ub)
        if isempty(lb); lb = -inf(N_u, 1); end
        if isempty(ub); ub =  inf(N_u, 1); end
        lb_z = [-inf(N_x, 1); lb];
        ub_z = [ inf(N_x, 1); ub];
    else
        lb_z = [];
        ub_z = [];
    end
    
elseif const_type == 1
    % Soft constraints boundary generation (includes eta bounds)
    if ~isempty(lb) || ~isempty(ub)
        if isempty(lb); lb = -inf(N_u, 1); end
        if isempty(ub); ub =  inf(N_u, 1); end
        % Enforce eta >= 0.0001
        lb_z = [-inf(N_x, 1); lb; ones(N_eta, 1) * 0.0001]; 
        ub_z = [ inf(N_x, 1); ub; inf(N_eta, 1)]; 
    else
        lb_z = [];
        ub_z = [];
    end
end

%% 6. Output Constraints (Inequalities)

if use_output_constraint
    if const_type == 0 
        % Hard Output Constraints: yhat(ii) <= y_max_inc
        C_til  = [zeros(H, n), kron(eye(H), C)];
        A_ineq = [C_til, zeros(H, N_u)];
        b_ineq = y_max_inc * ones(H, 1);
        
    elseif const_type == 1 
        % Soft Output Constraints: yhat(ii) - eta(ii) <= y_max_inc
        C_til  = [zeros(H, n), kron(eye(H), C)];
        A_ineq = [C_til, zeros(H, N_u), -eye(N_eta)]; 
        b_ineq = y_max_inc * ones(H, 1);
    end
else
    A_ineq = [];
    b_ineq = [];
end

%% 7. Solve Optimization Problem

% Configure Quadprog solver options
opts = optimoptions('quadprog', ...
    'Algorithm',           'interior-point-convex', ...
    'OptimalityTolerance', 1e-9,  ...
    'ConstraintTolerance', 1e-9,  ...
    'StepTolerance',       1e-12, ...
    'Display',             'off');

% Solve for the optimal decision vector
[U_opt, ~, exitflag] = quadprog(F, f, A_ineq, b_ineq, Aeq, beq, lb_z, ub_z, [], opts);

%% 8. Receding Horizon Step

% Extract and return only the first optimal control action
% U_opt contains states first (N_x), then control inputs (N_u)
u0 = U_opt(N_x + 1);

end