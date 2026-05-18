function [u0, exitflag] = mpc_solve(x0, H, R, A, B, C, lb, ub, y_max_inc)
% MPC_SOLVE  Solve the MPC regulator using quadprog (dense formulation).
%
%   [u0, exitflag] = mpc_solve(x0, H, R, A, B, C)
%   [u0, exitflag] = mpc_solve(x0, H, R, A, B, C, lb, ub)
%   [u0, exitflag] = mpc_solve(x0, H, R, A, B, C, lb, ub, y_max_inc)

%% ── Defaults ─────────────────────────────────────────────────────────────
if nargin < 7,  lb        = []; end
if nargin < 8,  ub        = []; end
if nargin < 9,  y_max_inc = []; end     % empty = no output constraint

use_output_constraint = ~isempty(y_max_inc);

%% ── Build prediction matrices ────────────────────────────────────────────
n = size(A, 1);

W  = zeros(H, H);
for ii = 1:H
    for jj = 1:ii
        W(ii, jj) = C * (A^(ii-jj)) * B;
    end
end

Pi = zeros(H, n);
for ii = 1:H
    Pi(ii, :) = C * (A^ii);
end

M = W'*W + R*eye(H);

%% ── quadprog cost matrices ───────────────────────────────────────────────
F = 2 * M;
f = 2 * W' * Pi * x0;

%% ── HARD CONSTRAINTS ──────────────────────────────────────────────────
if use_output_constraint
    A_ineq = W;
    b_ineq = y_max_inc * ones(H,1) - Pi * x0;
else
    A_ineq = [];
    b_ineq = [];
end


%% ── Solve ────────────────────────────────────────────────────────────────
opts = optimoptions('quadprog', 'Display', 'off');
[U_opt, ~, exitflag] = quadprog(F, f, A_ineq, b_ineq, [], [], lb, ub, [], opts);

%% ── Receding horizon: return only first control action ───────────────────
u0 = U_opt(1);

end