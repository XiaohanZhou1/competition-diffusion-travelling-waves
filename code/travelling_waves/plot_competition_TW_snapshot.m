function plot_competition_TW_snapshot()
% plot_competition_TW_snapshot_only
%
% Solve the two-species competition--diffusion system
%
%   u_t = u_xx + u(1-u)(u-a1) - uv/delta,
%   v_t = D v_xx + v(1-v)     - gamma*uv/delta,
%
% using a method-of-lines discretisation with ode15s.
%
% The final figure shows one snapshot at t_* = 40:
%
%   blue solid curve : u(x,t_*)
%   red solid curve  : v(x,t_*)
%   black dashed line: L(t_*) = argmax_x u(x,t_*)v(x,t_*)

    clearvars -except ans;
    clc;
    close all;

    %% ============================================================
    %  PARAMETERS
    % ============================================================

    % Model parameters
    delta = 1e-3;
    a1    = 0.30;
    D     = 5.0;
    gamma = 0.30;

    % Initial plateau heights
    d_u = 1.0;
    d_v = 1.0;

    % Domain / discretisation
    Xmax = 250;
    Nx   = 900;

    % Time
    tEnd  = 180;
    NtOut = 360;

    % Snapshot time
    t_snapshot = 40;

    % Initial separation
    gap = 60;

    % Smoothness of the smoothed Heaviside initial fronts
    wIC = 0.01;

    % Solver options
    RelTol  = 1e-5;
    AbsTol  = 1e-7;
    MaxStep = 0.5;

    % Interface tracking tolerance
    interaction_tol = 1e-6;

    %% ============================================================
    %  GRID AND INITIAL CONDITIONS
    % ============================================================

    x  = linspace(0, Xmax, Nx).';
    dx = x(2) - x(1);

    Lmat = neumann_laplacian_1d(Nx, dx);

    xmid = 0.5 * Xmax;

    % Initial front locations in the computational coordinate
    b_u = xmid - gap/2;
    b_v = xmid + gap/2;

    % Shifted plotting coordinate:
    % initial u-front at x = 0,
    % initial v-front at x = gap.
    x_plot = x - b_u;

    [u0, v0] = initial_condition_gap(x, d_u, d_v, wIC, b_u, b_v);
    y0 = [u0; v0];

    % Include t_snapshot exactly in the output times
    tspan = unique([linspace(0, tEnd, NtOut), t_snapshot]);

    %% ============================================================
    %  ODE15S SETUP
    % ============================================================

    rhs  = @(t,y) rhs_fun(t, y, Lmat, D, delta, gamma, a1);
    jac  = @(t,y) jacobian_fun(y, Lmat, D, delta, gamma, a1);
    Jpat = jacobian_pattern(Nx, Lmat);

    ode_opts = odeset( ...
        'RelTol', RelTol, ...
        'AbsTol', AbsTol, ...
        'MaxStep', MaxStep, ...
        'Jacobian', jac, ...
        'JPattern', Jpat, ...
        'NonNegative', 1:(2*Nx));

    %% ============================================================
    %  SOLVE
    % ============================================================

    fprintf('Solving PDE system with ode15s...\n');
    [t_sol, y_sol] = ode15s(rhs, tspan, y0, ode_opts);
    fprintf('Done.\n');

    u = y_sol(:, 1:Nx);
    v = y_sol(:, Nx+1:2*Nx);

    %% ============================================================
    %  TRACK INTERFACE L(t) = argmax_x u(x,t)v(x,t)
    % ============================================================

    Nt = numel(t_sol);

    uv_prod = u .* v;

    L_track = nan(Nt,1);
    max_uv  = nan(Nt,1);

    for k = 1:Nt
        [max_uv(k), idx_max] = max(uv_prod(k,:));

        if max_uv(k) > interaction_tol
            L_track(k) = x_plot(idx_max);
        end
    end

    %% ============================================================
    %  EXTRACT SNAPSHOT AT t_* = 40
    % ============================================================

    [~, k_snap] = min(abs(t_sol - t_snapshot));

    t_snap = t_sol(k_snap);
    u_snap = u(k_snap, :);
    v_snap = v(k_snap, :);
    L_snap = L_track(k_snap);

    fprintf('Snapshot time: t_* = %.4f\n', t_snap);
    fprintf('Interface location: L(t_*) = %.4f\n', L_snap);

    %% ============================================================
    %  CHOOSE PLOTTING WINDOW
    % ============================================================

    % Use a window centred around the interaction region,
    % but also wide enough to show both profiles clearly.
    if isfinite(L_snap)
        x_left  = L_snap - 70;
        x_right = L_snap + 90;
    else
        x_left  = -80;
        x_right = 80;
    end

    x_left  = max(x_left,  min(x_plot));
    x_right = min(x_right, max(x_plot));

    idx_xwin = (x_plot >= x_left) & (x_plot <= x_right);

    x_fig = x_plot(idx_xwin);
    u_fig = u_snap(idx_xwin);
    v_fig = v_snap(idx_xwin);

    %% ============================================================
    %  PLOT SNAPSHOT PROFILES
    % ============================================================

    figure('Color','w', 'Position', [100, 100, 650, 450]);
    hold on;

    plot(x_fig, u_fig, '-', ...
        'Color', [0 0.25 0.9], ...
        'LineWidth', 2.3, ...
        'DisplayName', '$u(x,t_*)$');

    plot(x_fig, v_fig, '-', ...
        'Color', [0.9 0 0], ...
        'LineWidth', 2.3, ...
        'DisplayName', '$v(x,t_*)$');

    if isfinite(L_snap)
        plot([L_snap, L_snap], [0, 1.05], 'k--', ...
            'LineWidth', 2.0, ...
            'DisplayName', '$L(t_*)$');
    end

    % Text box showing the snapshot time
    text(0.05, 0.92, '$t_*=40$', ...
        'Units', 'normalized', ...
        'Interpreter', 'latex', ...
        'FontSize', 15, ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.4 0.4 0.4], ...
        'Margin', 6);

    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Population density', 'Interpreter', 'latex', 'FontSize', 18);

    xlim([x_left, x_right]);
    ylim([0, 1.05]);

    set(gca, ...
        'FontSize', 14, ...
        'LineWidth', 1.1, ...
        'TickLabelInterpreter', 'latex');

    box on;
    grid on;

    legend('Interpreter', 'latex', ...
        'FontSize', 12, ...
        'Location', 'northeast');

    %% ============================================================
    %  OPTIONAL EXPORT
    % ============================================================

    % exportgraphics(gcf, 'competition_TW_snapshot_t40.png', 'Resolution', 300);
    % exportgraphics(gcf, 'competition_TW_snapshot_t40.pdf', 'ContentType', 'vector');

end

%% ========================================================================
%  RHS
% ========================================================================
function dydt = rhs_fun(~, y, Lmat, D, delta, gamma, a1)

    Nx = numel(y)/2;

    u = y(1:Nx);
    v = y(Nx+1:end);

    R = u .* (1-u) .* (u-a1);
    S = v .* (1-v);

    dudt = Lmat*u + R - (u.*v)/delta;
    dvdt = D*(Lmat*v) + S - gamma*(u.*v)/delta;

    dydt = [dudt; dvdt];
end

%% ========================================================================
%  SPARSE JACOBIAN
% ========================================================================
function J = jacobian_fun(y, Lmat, D, delta, gamma, a1)

    Nx = numel(y)/2;

    u = y(1:Nx);
    v = y(Nx+1:end);

    Rp = -3*u.^2 + 2*(1+a1)*u - a1;
    Sp = 1 - 2*v;

    J11 = Lmat + spdiags(Rp - v/delta, 0, Nx, Nx);
    J12 = spdiags(-u/delta, 0, Nx, Nx);

    J21 = spdiags(-gamma*v/delta, 0, Nx, Nx);
    J22 = D*Lmat + spdiags(Sp - gamma*u/delta, 0, Nx, Nx);

    J = [J11, J12;
         J21, J22];
end

%% ========================================================================
%  JACOBIAN SPARSITY PATTERN
% ========================================================================
function Jpat = jacobian_pattern(Nx, Lmat)

    P = spones(Lmat);
    I = speye(Nx);

    Jpat = [P + I, I;
            I,     P + I];
end

%% ========================================================================
%  NEUMANN LAPLACIAN
% ========================================================================
function L = neumann_laplacian_1d(Nx, dx)

    e = ones(Nx,1);
    L = spdiags([e -2*e e], -1:1, Nx, Nx);

    % Homogeneous Neumann BC via ghost-point discretisation
    L(1,2)     = 2;
    L(Nx,Nx-1) = 2;

    L = L / dx^2;
end

%% ========================================================================
%  INITIAL CONDITIONS
% ========================================================================
function [u0, v0] = initial_condition_gap(x, d_u, d_v, w, b_u, b_v)

    H1 = 0.5 * (1 + tanh((x - b_u)/w));
    H2 = 0.5 * (1 + tanh((x - b_v)/w));

    u0 = d_u * (1 - H1);
    v0 = d_v * H2;
end