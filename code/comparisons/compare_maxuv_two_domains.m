function compare_maxuv_two_domains()
% compare_maxuv_two_domains_semilogy
%
% Compare M(t) = max_x u(x,t)v(x,t) for two computational domains.
%
% The purpose is to test whether increasing the x-domain delays the
% late-time loss of the stable overlap regime, while leaving the plateau
% value of max_x uv essentially unchanged.
%
% Model:
%   u_t = u_xx + u(1-u)(u-a1) - uv/delta,
%   v_t = D v_xx + v(1-v)     - gamma uv/delta.
%
% Numerical method:
%   Method of lines + ode15s with sparse Jacobian.

    clearvars -except ans;
    clc;
    close all;

    %% ============================================================
    %  MODEL PARAMETERS
    % ============================================================

    delta = 1e-3;
    a1    = 0.20;
    D     = 5.0;
    gamma = 0.30;

    % Initial plateau heights
    d_u = 1.0;
    d_v = 1.0;

    % Initial separation
    gap = 60;

    % Smoothness of initial fronts
    wIC = 0.01;

    %% ============================================================
    %  DOMAIN AND TIME PARAMETERS
    % ============================================================

    % Compare a shorter and a longer computational domain
    Xmax_list = [250, 400];

    % Use approximately the same dx as in the Xmax = 250, Nx = 900 case
    dx_ref = 250 / (900 - 1);

    % Time interval
    tEnd  = 500;
    NtOut = 900;
    tspan = linspace(0, tEnd, NtOut);

    %% ============================================================
    %  SOLVER OPTIONS
    % ============================================================

    RelTol  = 1e-5;
    AbsTol  = 1e-7;
    MaxStep = 0.5;

    %% ============================================================
    %  DIAGNOSTIC OPTIONS
    % ============================================================

    % Plotting floor for semilogy only.
    % This does not change the computed maxuv values.
    plot_floor = 1e-25;

    % Threshold used to identify effective loss of overlap.
    overlap_tol = 1e-13;

    % Stable time window for comparing plateau values.
    % Adjust this if needed for your parameter regime.
    stable_window = [30, 70];

    %% ============================================================
    %  RUN SIMULATIONS
    % ============================================================

    nDom = numel(Xmax_list);
    results = struct();

    for j = 1:nDom

        Xmax = Xmax_list(j);

        % Keep dx approximately fixed when Xmax changes
        Nx = round(Xmax / dx_ref) + 1;

        fprintf('\n============================================================\n');
        fprintf('Running domain %d/%d: Xmax = %.1f, Nx = %d\n', ...
            j, nDom, Xmax, Nx);
        fprintf('============================================================\n');

        [t_sol, x_plot, u, v, L_track, maxuv] = run_one_domain( ...
            Xmax, Nx, tspan, ...
            delta, a1, D, gamma, d_u, d_v, gap, wIC, ...
            RelTol, AbsTol, MaxStep);

        results(j).Xmax    = Xmax;
        results(j).Nx      = Nx;
        results(j).t       = t_sol;
        results(j).x_plot  = x_plot;
        results(j).u       = u;
        results(j).v       = v;
        results(j).L_track = L_track;
        results(j).maxuv   = maxuv;

        % Find first time after contact/plateau when maxuv falls below threshold
        idx_after_stable = t_sol > stable_window(2);
        idx_loss = find(idx_after_stable & (maxuv < overlap_tol), 1, 'first');

        if ~isempty(idx_loss)
            results(j).loss_time = t_sol(idx_loss);
            fprintf('Loss of overlap below %.1e occurs at t ~= %.4f\n', ...
                overlap_tol, results(j).loss_time);
        else
            results(j).loss_time = NaN;
            fprintf('No loss of overlap below %.1e detected up to tEnd.\n', overlap_tol);
        end

        % Stable-window statistics
        idx_stable = (t_sol >= stable_window(1)) & (t_sol <= stable_window(2));

        mean_maxuv = mean(maxuv(idx_stable), 'omitnan');
        std_maxuv  = std(maxuv(idx_stable),  'omitnan');

        results(j).mean_maxuv = mean_maxuv;
        results(j).std_maxuv  = std_maxuv;

        fprintf('Stable window t in [%.1f, %.1f]: mean maxuv = %.6e, std = %.6e\n', ...
            stable_window(1), stable_window(2), mean_maxuv, std_maxuv);
    end

    %% ============================================================
    %  PLOT: MAXIMUM OVERLAP ON SEMILOGY SCALE
    % ============================================================

    figure('Color','w', 'Position', [100, 100, 760, 500]);
    hold on;

    line_styles = {'-', '--', '-.', ':'};

    for j = 1:nDom

        t = results(j).t;
        maxuv = results(j).maxuv;

        % Plotting only: replace zeros and tiny values by a floor
        maxuv_plot = max(maxuv, plot_floor);

        semilogy(t, maxuv_plot, ...
            'LineStyle', line_styles{j}, ...
            'LineWidth', 2.2, ...
            'DisplayName', sprintf('$X_{\\max}=%.0f$', results(j).Xmax));
    end

    yline(overlap_tol, 'k--', ...
        'LineWidth', 1.3, ...
        'DisplayName', sprintf('$%.0e$', overlap_tol));

    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$\max_x u(x,t)v(x,t)$', ...
        'Interpreter', 'latex', 'FontSize', 18);

    title('Maximum overlap', ...
        'Interpreter', 'latex', 'FontSize', 16);

    legend('Interpreter', 'latex', ...
        'FontSize', 12, ...
        'Location', 'best');

    set(gca, ...
        'FontSize', 14, ...
        'LineWidth', 1.1, ...
        'TickLabelInterpreter', 'latex', ...
        'YScale', 'log');

    ylim([plot_floor, 1]);
    xlim([0, tEnd]);

    box on;
    grid on;

    %% ============================================================
    %  OPTIONAL: PRINT RELATIVE DIFFERENCE IN PLATEAU VALUES
    % ============================================================

    if nDom == 2
        rel_diff = abs(results(1).mean_maxuv - results(2).mean_maxuv) ...
                   / max(results(1).mean_maxuv, results(2).mean_maxuv);

        fprintf('\nRelative difference in stable mean maxuv: %.4e\n', rel_diff);
    end

    %% ============================================================
    %  OPTIONAL EXPORT
    % ============================================================

    % exportgraphics(gcf, 'compare_maxuv_two_domains_semilogy.png', 'Resolution', 300);
    % exportgraphics(gcf, 'compare_maxuv_two_domains_semilogy.pdf', 'ContentType', 'vector');

end

%% ========================================================================
%  RUN ONE DOMAIN
% ========================================================================
function [t_sol, x_plot, u, v, L_track, maxuv] = run_one_domain( ...
    Xmax, Nx, tspan, ...
    delta, a1, D, gamma, d_u, d_v, gap, wIC, ...
    RelTol, AbsTol, MaxStep)

    %% Grid

    x  = linspace(0, Xmax, Nx).';
    dx = x(2) - x(1);

    Lmat = neumann_laplacian_1d(Nx, dx);

    xmid = 0.5 * Xmax;

    % Initial front locations in computational coordinate
    b_u = xmid - gap/2;
    b_v = xmid + gap/2;

    % Shifted coordinate:
    % initial u-front at x = 0,
    % initial v-front at x = gap.
    x_plot = x - b_u;

    [u0, v0] = initial_condition_gap(x, d_u, d_v, wIC, b_u, b_v);
    y0 = [u0; v0];

    %% ode15s setup

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

    %% Solve

    fprintf('Solving with ode15s...\n');
    [t_sol, y_sol] = ode15s(rhs, tspan, y0, ode_opts);
    fprintf('Done.\n');

    u = y_sol(:, 1:Nx);
    v = y_sol(:, Nx+1:2*Nx);

    %% Diagnostics

    uv_prod = u .* v;

    % Raw maximum overlap. Do not threshold here.
    maxuv = max(uv_prod, [], 2);

    Nt = numel(t_sol);
    L_track = nan(Nt,1);

    interaction_tol = 1e-14;

    for k = 1:Nt
        [mval, idx_max] = max(uv_prod(k,:));

        if mval > interaction_tol
            L_track(k) = x_plot(idx_max);
        end
    end

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