function plot_competition_TW_profiles()
% plot_competition_TW_topdown_single
%
% Solve the two-species competition--diffusion system
%
%   u_t = u_xx + u(1-u)(u-a1) - uv/delta,
%   v_t = D v_xx + v(1-v)     - gamma*uv/delta,
%
% using a method-of-lines discretisation with ode15s.
%
% Plot one top-down surface plot:
%   - background colour map: u(x,t)
%   - solid contours        : selected levels of v(x,t)
%   - light dashed curve    : L(t) = argmax_x u(x,t)v(x,t)
%
% Axes:
%   x-axis : space x
%   y-axis : time t

    clearvars -except ans;
    clc;
    close all;

    %% ============================================================
    %  PARAMETERS
    % ============================================================

    % Model parameters
    delta = 1e-3;
    a1    = 0.20;
    D     = 5.0;
    gamma = 0.30;

    % Initial plateau heights
    d_u = 0.5;
    d_v = 0.5;

    % Domain / discretisation
    Xmax = 250;
    Nx   = 900;

    % Time
    tEnd  = 180;
    NtOut = 360;

    % Initial separation
    gap = 60;

    % Smoothness of the Heaviside approximation
    wIC = 0.01;

    % Solver options
    RelTol  = 1e-5;
    AbsTol  = 1e-7;
    MaxStep = 0.5;

    % Interface tracking tolerance:
    % only plot L(t) when max_x u(x,t)v(x,t) exceeds this threshold
    interaction_tol = 1e-6;

    % Time window shown in the final figure
    t_plot_max = 90;
    % If you want the whole simulation, set t_plot_max = tEnd;

    % Contour levels for v(x,t)
    v_levels = [0.1 0.3 0.5 0.7 0.9];

    %% ============================================================
    %  GRID AND INITIAL CONDITIONS
    % ============================================================

    x  = linspace(0, Xmax, Nx).';
    dx = x(2) - x(1);

    Lmat = neumann_laplacian_1d(Nx, dx);

    xmid = 0.5 * Xmax;

    % In the computational coordinate:
    b_u = xmid - gap/2;
    b_v = xmid + gap/2;

    % Shifted plotting coordinate:
    %   initial u-front at x = 0
    %   initial v-front at x = gap
    x_plot = x - b_u;

    [u0, v0] = initial_condition_gap(x, d_u, d_v, wIC, b_u, b_v);
    y0 = [u0; v0];

    tspan = linspace(0, tEnd, NtOut);

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

    k_first_interaction = find(isfinite(L_track), 1, 'first');
    if ~isempty(k_first_interaction)
        fprintf('First plotted interaction time: t ~= %.4f\n', t_sol(k_first_interaction));
        fprintf('Initial plotted interface location: L(t) ~= %.4f\n', L_track(k_first_interaction));
    else
        warning('No interaction region detected with the current interaction_tol.');
    end

    %% ============================================================
    %  CHOOSE SPATIAL WINDOW
    % ============================================================

    if any(isfinite(L_track))
        L_min = min(L_track(isfinite(L_track)));
        L_max = max(L_track(isfinite(L_track)));

        x_left  = L_min - 45;
        x_right = L_max + 45;
    else
        x_left  = -60;
        x_right = gap + 120;
    end

    % Also keep the initial fronts visible
    x_left  = min(x_left, -60);
    x_right = max(x_right, gap + 60);

    % Restrict to computational domain
    x_left  = max(x_left,  min(x_plot));
    x_right = min(x_right, max(x_plot));

    idx_xwin = (x_plot >= x_left) & (x_plot <= x_right);

    %% ============================================================
    %  CHOOSE TIME WINDOW
    % ============================================================

    idx_twin = (t_sol <= t_plot_max);

    t_fig = t_sol(idx_twin);
    u_fig = u(idx_twin, idx_xwin);
    v_fig = v(idx_twin, idx_xwin);
    L_fig = L_track(idx_twin);

    x_fig = x_plot(idx_xwin);

    density_max = max([max(u_fig, [], 'all'), max(v_fig, [], 'all'), 1.0]);

    %% ============================================================
    %  PLOT: u(x,t) BACKGROUND + v(x,t) CONTOURS + L(t)
    % ============================================================

    figure('Color','w', 'Position', [100, 100, 760, 500]);
    hold on;

    % Background colour map: u(x,t)
    surf(x_fig, t_fig, u_fig, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');
    view(2);
    shading interp;
    colormap(parula);
    caxis([0, density_max]);

    % Contours of v(x,t): SOLID LINES
    [~, hc] = contour(x_fig, t_fig, v_fig, v_levels, ...
        'LineColor', 'w', ...
        'LineStyle', '-', ...
        'LineWidth', 1.1);

    % Interface L(t): LIGHT DASHED LINE
    validL = isfinite(L_fig);
    plot(L_fig(validL), t_fig(validL), ...
        '--', ...
        'Color', [0.78 0.78 0.78], ...
        'LineWidth', 2.2);

    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);

    title('$u(x,t)$ with contours of $v(x,t)$ and interface $L(t)$', ...
        'Interpreter', 'latex', 'FontSize', 17);

    cb = colorbar;
    cb.Label.String = '$u(x,t)$';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';

    xlim([x_left, x_right]);
    ylim([min(t_fig), max(t_fig)]);

    set(gca, ...
        'FontSize', 14, ...
        'LineWidth', 1.1, ...
        'TickLabelInterpreter', 'latex', ...
        'YDir', 'normal');

    box on;

    % Optional legend
    h1 = plot(nan, nan, '-', ...
        'Color', 'w', ...
        'LineWidth', 1.2);

    h2 = plot(nan, nan, '--', ...
        'Color', [0.78 0.78 0.78], ...
        'LineWidth', 2.2);

    legend([h1, h2], ...
        {'Contours of $v(x,t)$', '$L(t)=\arg\max_x u(x,t)v(x,t)$'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 12, ...
        'Location', 'northeast');

    %% ============================================================
    %  OPTIONAL: EXPORT FIGURE
    % ============================================================

    % Uncomment if needed:
    % exportgraphics(gcf, 'competition_TW_topdown_single.png', 'Resolution', 300);
    % exportgraphics(gcf, 'competition_TW_topdown_single.pdf', 'ContentType', 'vector');

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