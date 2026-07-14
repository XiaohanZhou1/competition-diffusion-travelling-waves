function plot_competition_TW_RGB_composite()
% plot_competition_TW_RGB_composite
%
% Solve the two-species competition--diffusion system
%
%   u_t = u_xx + u(1-u)(u-a1) - uv/delta,
%   v_t = D v_xx + v(1-v)     - gamma*uv/delta,
%
% using a method-of-lines discretisation with ode15s.
%
% The final figure is a top-down RGB composite plot:
%
%   blue  : regions dominated by u(x,t),
%   red   : regions dominated by v(x,t),
%   white : regions where both densities are small.
%
% The dashed curve denotes the computed interface location
%
%   L(t) = argmax_x u(x,t)v(x,t).

    clearvars -except ans;
    clc;
    close all;

    %% ============================================================
    %  PARAMETERS
    % ============================================================

    % Model parameters
    delta = 1e-4;
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

    % Initial separation
    gap = 60;

    % Smoothness of the smoothed Heaviside initial fronts
    wIC = 0.01;

    % Solver options
    RelTol  = 1e-5;
    AbsTol  = 1e-7;
    MaxStep = 0.5;

    % Interface tracking tolerance:
    % L(t) is plotted only when max_x u(x,t)v(x,t) exceeds this threshold.
    interaction_tol = 1e-6;

    % Time window shown in the final figure.
    % Increase to tEnd if you want to show the whole simulation.
    t_plot_min = 0;
    t_plot_max = 90;

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
    %   initial u-front is at x = 0,
    %   initial v-front is at x = gap.
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
    %  CHOOSE TIME WINDOW
    % ============================================================

    idx_twin = (t_sol >= t_plot_min) & (t_sol <= t_plot_max);

    t_fig = t_sol(idx_twin);
    u_twin = u(idx_twin, :);
    v_twin = v(idx_twin, :);
    L_fig = L_track(idx_twin);

    %% ============================================================
    %  CHOOSE SPATIAL WINDOW
    % ============================================================

    if any(isfinite(L_fig))
        L_min = min(L_fig(isfinite(L_fig)));
        L_max = max(L_fig(isfinite(L_fig)));

        x_left  = L_min - 45;
        x_right = L_max + 45;
    else
        x_left  = -60;
        x_right = gap + 80;
    end

    % Keep the initial gap and initial v-front visible
    x_left  = min(x_left, -70);
    x_right = max(x_right, gap + 70);

    % Restrict to computational domain
    x_left  = max(x_left,  min(x_plot));
    x_right = min(x_right, max(x_plot));

    idx_xwin = (x_plot >= x_left) & (x_plot <= x_right);

    x_fig = x_plot(idx_xwin);
    U = u_twin(:, idx_xwin);
    V = v_twin(:, idx_xwin);



    %% ============================================================
    %  RGB COMPOSITE COLOUR FIELD WITH LOW-DENSITY CUTOFF
    % ============================================================
    %
    % The colour field is constructed from thresholded densities.
    % Values below eta_colour are treated as zero so that very small
    % diffusive tails appear as low-density/white regions.
    
    eta_colour = 0;
    
    Uplot = U;
    Vplot = V;
    
    % Remove very small diffusive tails from the colour representation
    Uplot(Uplot < eta_colour) = 0;
    Vplot(Vplot < eta_colour) = 0;
    
    % Keep values in [0,1]
    Uplot = min(max(Uplot, 0), 1);
    Vplot = min(max(Vplot, 0), 1);
    
    ColorRedBlue = zeros(size(Uplot,1), size(Uplot,2), 3);
    
    ColorRedBlue(:,:,1) = 1 - Uplot;
    ColorRedBlue(:,:,2) = 1 - (Uplot + Vplot);
    ColorRedBlue(:,:,3) = 1 - Vplot;
    
    ColorRedBlue = max(0, min(1, ColorRedBlue));

    %% ============================================================
    %  PLOT RGB COMPOSITE FIELD + INTERFACE L(t)
    % ============================================================

    [Xgrid, Tgrid] = meshgrid(x_fig, t_fig);
    Zgrid = zeros(size(U));

    figure('Color','w', 'Position', [100, 100, 760, 500]);
    hold on;

    surf(Xgrid, Tgrid, Zgrid, ColorRedBlue, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'flat');

    view(2);
    axis tight;

    % Computed interface location L(t)
    validL = isfinite(L_fig);
    plot(L_fig(validL), t_fig(validL), '--', ...
        'Color', [0.15 0.15 0.15], ...
        'LineWidth', 2.2);

    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);

    xlim([x_left, x_right]);
    ylim([min(t_fig), max(t_fig)]);

    set(gca, ...
        'FontSize', 14, ...
        'LineWidth', 1.1, ...
        'TickLabelInterpreter', 'latex', ...
        'YDir', 'normal');

    box on;

    %% ============================================================
    %  LEGEND
    % ============================================================

    h_u = plot(nan, nan, 's', ...
        'MarkerFaceColor', [0 0 1], ...
        'MarkerEdgeColor', [0 0 1], ...
        'MarkerSize', 8);

    h_v = plot(nan, nan, 's', ...
        'MarkerFaceColor', [1 0 0], ...
        'MarkerEdgeColor', [1 0 0], ...
        'MarkerSize', 8);

    h_n = plot(nan, nan, 's', ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerSize', 8);

    h_L = plot(nan, nan, '--', ...
        'Color', [0.15 0.15 0.15], ...
        'LineWidth', 2.2);

    legend([h_u, h_v, h_n, h_L], ...
        {'$u$-dominated', '$v$-dominated', 'low density', '$L(t)$'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 12, ...
        'Location', 'northeast');

    %% ============================================================
    %  OPTIONAL: EXPORT FIGURE
    % ============================================================

    % Uncomment if needed:
    %
    % exportgraphics(gcf, 'competition_TW_RGB_composite.png', 'Resolution', 300);
    % exportgraphics(gcf, 'competition_TW_RGB_composite.pdf', 'ContentType', 'vector');

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