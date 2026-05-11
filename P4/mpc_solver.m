function u0 = mpc_solve(x0, H, R, A, B, C)
% MPC_SOLVE  Solve the unconstrained MPC regulator using quadprog.
%
%   u0 = mpc_solve(x0, H, R, A, B, C)
%
%   Solves the finite-horizon quadratic MPC problem (dense formulation):
%
%       min_U  J = Y'Y + R*U'U
%
%   where Y = W*U + Pi*x0, which expands to:
%
%       min_U  U'*F/2*U + f'*U
%
%   with:
%       F = 2*M = 2*(W'W + R*I)     [H x H]
%       f = 2*W'*Pi*x0              [H x 1]
%
%   Only the first element of U* is returned (receding horizon).
%
%   Inputs
%     x0  – n×1 current state (incremental: Dx at current time)
%     H   – prediction horizon (positive integer)
%     R   – positive scalar control weight
%     A   – n×n state matrix  (from identified model)
%     B   – n×1 input matrix  (scalar input assumed)
%     C   – 1×n output matrix (scalar output assumed)
%
%   Output
%     u0  – scalar, first optimal control action to apply to the plant
%
%   Note: this function will be extended in Q3 (control constraints),
%   Q4 (tracking via change of variables) and Q5 (soft output constraint).

%% ── Build prediction matrices (same as compute_KRH) ─────────────────────
n = size(A, 1);

% W  (H×H) lower-triangular: W(i,j) = C * A^(i-j) * B  for i >= j
W = zeros(H, H);
for ii = 1:H
    for jj = 1:ii
        W(ii, jj) = C * (A^(ii-jj)) * B;
    end
end

% Pi  (H×n): Pi(i,:) = C * A^i
Pi = zeros(H, n);
for ii = 1:H
    Pi(ii, :) = C * (A^ii);
end

% Hessian of the cost w.r.t. U
M = W'*W + R*eye(H);

%% ── quadprog cost matrices (dense formulation) ───────────────────────────
% J = U'*M*U + 2*x0'*Pi'*W*U + const
% quadprog minimises  1/2 * z'*F*z + f'*z
% so F = 2*M,  f = 2*W'*Pi*x0
F = 2 * M;
f = 2 * W' * Pi * x0;

%% ── Solve with quadprog ──────────────────────────────────────────────────
opts = optimoptions('quadprog', 'Display', 'off');
U_opt = quadprog(F, f, [], [], [], [], [], [], [], opts);

%% ── Receding horizon: apply only first control action ────────────────────
u0 = U_opt(1);

end