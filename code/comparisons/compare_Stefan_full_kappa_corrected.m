function compare_Stefan_full_kappa_corrected()
% compare_Stefan_full_kappa_corrected
%
% Compare the full two-species PDE with the reduced Stefan-type travelling
% wave model using the corrected effective Stefan parameter.
%
% Convention used here:
%   The reduced Stefan-type curve is the phase-plane curve
%       kappa = -c/W(0;c).
%   The full two-species speed is c = dL/dt, and its horizontal coordinate
%   is computed from
%       kappa_corrected = 2*gamma*c^2/(2*c^2 + D).

    clearvars -except ans;
    clc;
    close all;

    %% ============================================================
    %  USER PARAMETERS
    % ============================================================
    D_fixed = 0.01;
    a       = 0.30;
    delta   = 1e-4;

    nFull = 31;
    gamma_list = linspace(0.5, 25, nFull);
    plot_positive_speed_only = true;
    min_R2_for_plot = -Inf;

    %% ============================================================
    %  REDUCED STEFAN-TYPE MODEL
    % ============================================================
    stefan = compute_stefan_curve(a, D_fixed);

    %% ============================================================
    %  FULL TWO-SPECIES PDE
    % ============================================================
    opts = default_full_pde_options();
    opts.delta = delta;
    opts.a1 = a;

    rows = nan(numel(gamma_list), 7);

    fprintf('Full PDE vs reduced Stefan comparison\n');
    fprintf('  D = %.6g\n', D_fixed);
    fprintf('  a = %.6g\n', a);
    fprintf('  x-axis: kappa_corrected = 2*gamma*c^2/(2*c^2 + D)\n');
    fprintf('  plotted two-species speed is c = dL/dt\n');
    fprintf('  full PDE points = %d\n\n', numel(gamma_list));

    for i = 1:numel(gamma_list)
        gamma = gamma_list(i);

        fprintf('Full PDE case %d/%d: gamma = %.8f\n', ...
            i, numel(gamma_list), gamma);

        out = run_one_full_pde_case(D_fixed, gamma, opts);

        c_full = out.c;
        kappa_full = corrected_kappa(gamma, c_full, D_fixed);
        rows(i,:) = [D_fixed, gamma, kappa_full, c_full, out.slope, out.R2, out.nFit];

        fprintf('  c = dL/dt = %+.8e, kappa_corrected = %.8f, R2 = %.4f, nFit = %d, method = %s\n', ...
            c_full, kappa_full, out.R2, out.nFit, out.fit_method);
    end

    Tfull = array2table(rows, 'VariableNames', ...
        {'D','gamma','kappa_corrected','c','slopeL','R2','nFit'});
    Tfull = sortrows(Tfull, 'kappa_corrected');

    Tplot = Tfull(isfinite(Tfull.kappa_corrected) & isfinite(Tfull.c), :);
    if plot_positive_speed_only
        Tplot = Tplot(Tplot.c > 0, :);
    end
    if isfinite(min_R2_for_plot)
        Tplot = Tplot(Tplot.R2 >= min_R2_for_plot, :);
    end

    %% ============================================================
    %  PLOT: c versus corrected kappa
    % ============================================================
    fig = figure('Color', 'w', 'Position', [120 120 720 520]);
    ax = axes(fig);
    hold(ax, 'on');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.65 0.65 0.65], 'MinorGridColor', [0.82 0.82 0.82]);

    hFull = plot(ax, Tplot.kappa_corrected, Tplot.c, '-', ...
        'Color', [0.85 0.10 0.05], ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Two-species model');

    kappa_plot_max = 1.08 * max([Tplot.kappa_corrected; stefan.kappa(:)]);
    stefan_plot = stefan.kappa <= kappa_plot_max;

    hStefan = plot(ax, stefan.kappa(stefan_plot), stefan.c(stefan_plot), '-', ...
        'Color', [0.0 0.25 1.0], ...
        'LineWidth', 2.0, ...
        'DisplayName', 'Reduced Stefan-type model');

    xlabel(ax, '$\kappa_{\mathrm{corrected}}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(ax, '$c$', 'Interpreter', 'latex', 'FontSize', 20);
    lgd = legend(ax, [hFull, hStefan], 'Location', 'best', 'Interpreter', 'latex');
    lgd.FontSize = 12;
    lgd.Color = 'w';
    lgd.TextColor = 'k';
    lgd.EdgeColor = [0.2 0.2 0.2];
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

    xlim(ax, [0, kappa_plot_max]);
    ylim(ax, [min([-0.05; Tplot.c(:)]), 1.05*max([stefan.c(:); Tplot.c(:); 0.05])]);
end

%% ========================================================================
%  FULL PDE OPTIONS
% ========================================================================
function opts = default_full_pde_options()
    opts = struct();
    opts.delta = 1e-4;
    opts.a1    = 0.30;
    opts.A     = 0.5;

    opts.Xmax = 120;
    opts.Nx   = 6000;
    opts.tEnd  = 160;
    opts.NtOut = 801;

    opts.gap = 60;
    opts.wIC = 0.01;

    opts.overlap_tol = 1e-6;
    opts.t_fit_start = 105;
    opts.fit_t_min = 105;
    opts.fit_t_max = 135;
    opts.fit_last_fraction = 0.40;

    opts.ic_type = 'gap';
    opts.interface_for_speed = 'uvmax';
    opts.local_centroid_half_window = 14;

    opts.RelTol  = 1e-5;
    opts.AbsTol  = 1e-7;
    opts.MaxStep = 0.25;

    opts.use_median_speed_fallback = false;
    opts.R2_fallback_threshold = 0.95;
end

%% ========================================================================
%  SOLVE ONE FULL PDE CASE
% ========================================================================
function out = run_one_full_pde_case(D, gamma, opts)
    delta = pick(opts, 'delta', 1e-3);
    a1    = pick(opts, 'a1', 0.30);
    A     = pick(opts, 'A', 1.0);

    Xmax = pick(opts, 'Xmax', 300);
    Nx   = pick(opts, 'Nx', 650);

    tEnd  = pick(opts, 'tEnd', 220);
    NtOut = pick(opts, 'NtOut', 360);

    gap = pick(opts, 'gap', 80);
    wIC = pick(opts, 'wIC', 1.0);
    ic_type = pick(opts, 'ic_type', 'gap');

    overlap_tol = pick(opts, 'overlap_tol', 1e-9);
    t_fit_start = pick(opts, 't_fit_start', 20);
    fit_t_min = pick(opts, 'fit_t_min', NaN);
    fit_t_max = pick(opts, 'fit_t_max', NaN);
    fit_last_fraction = pick(opts, 'fit_last_fraction', 0.40);

    interface_for_speed = pick(opts, 'interface_for_speed', 'localcentroid');
    local_win = pick(opts, 'local_centroid_half_window', 14);

    RelTol  = pick(opts, 'RelTol', 1e-5);
    AbsTol  = pick(opts, 'AbsTol', 1e-7);
    MaxStep = pick(opts, 'MaxStep', 0.8);

    use_median_speed_fallback = pick(opts, 'use_median_speed_fallback', true);
    R2_fallback_threshold     = pick(opts, 'R2_fallback_threshold', 0.95);

    x  = linspace(0, Xmax, Nx).';
    dx = x(2) - x(1);
    Lmat = neumann_laplacian_1d(Nx, dx);

    xmid    = 0.5 * Xmax;
    x_uedge = xmid - gap/2;
    x_vedge = xmid + gap/2;

    switch lower(ic_type)
        case 'gap'
            [u0, v0] = initial_condition_gap(x, A, wIC, x_uedge, x_vedge);
        case 'stefan_small_branch'
            [u0, v0] = initial_condition_stefan_small_branch(x, A, wIC, xmid, a1, D);
        otherwise
            error('Unknown ic_type = %s', ic_type);
    end
    y0 = [u0; v0];

    tspan = linspace(0, tEnd, NtOut);
    rhs = @(t,y) rhs_fun(t, y, Lmat, D, delta, gamma, a1);
    jac = @(t,y) jacobian_fun(y, Lmat, D, delta, gamma, a1);
    Jpat = jacobian_pattern(Nx, Lmat);

    ode_opts = odeset( ...
        'RelTol', RelTol, ...
        'AbsTol', AbsTol, ...
        'MaxStep', MaxStep, ...
        'Jacobian', jac, ...
        'JPattern', Jpat, ...
        'NonNegative', 1:(2*Nx));

    [t_sol, y_sol] = ode15s(rhs, tspan, y0, ode_opts);

    Nt = numel(t_sol);
    u = y_sol(:, 1:Nx);
    v = y_sol(:, Nx+1:2*Nx);

    overlap_max  = nan(Nt,1);
    L_uvmax      = nan(Nt,1);
    L_uvlocal    = nan(Nt,1);
    L_uvcentroid = nan(Nt,1);

    for k = 1:Nt
        uk = u(k,:).';
        vk = v(k,:).';

        wk = uk .* vk;
        [wmax, idxMax] = max(wk);
        overlap_max(k) = wmax;

        if wmax > overlap_tol
            L_uvmax(k) = subgrid_peak_location(x, wk, idxMax);

            i1 = max(1, idxMax - local_win);
            i2 = min(Nx, idxMax + local_win);
            x_loc = x(i1:i2);
            w_loc = wk(i1:i2);
            wloc_int = trapz(x_loc, w_loc);

            if isfinite(wloc_int) && wloc_int > overlap_tol
                L_uvlocal(k) = trapz(x_loc, x_loc .* w_loc) / wloc_int;
            end
        end

        w_int = trapz(x, wk);
        if isfinite(w_int) && w_int > overlap_tol
            L_uvcentroid(k) = trapz(x, x .* wk) / w_int;
        end
    end

    idxInteract = find(overlap_max > overlap_tol, 1, 'first');
    if isempty(idxInteract)
        fit_start_eff = t_fit_start;
    else
        fit_start_eff = max(t_fit_start, t_sol(idxInteract));
    end

    switch lower(interface_for_speed)
        case 'uvmax'
            Lfit = L_uvmax;
        case 'localcentroid'
            Lfit = L_uvlocal;
        case 'centroid'
            Lfit = L_uvcentroid;
        case 'reaction_peak'
            Lfit = L_uvmax;
        otherwise
            error('Unknown interface_for_speed = %s', interface_for_speed);
    end

    if isfinite(fit_t_min) && isfinite(fit_t_max)
        fit_idx = find(t_sol >= fit_t_min & t_sol <= fit_t_max & isfinite(Lfit));
        fit = fit_one_window(t_sol, Lfit, fit_idx);
        fit_method = 'fixed-window-linear';
    else
        fit = fit_interface_speed_active_window(t_sol, Lfit, overlap_max, overlap_tol, ...
                                                fit_start_eff, fit_last_fraction);
        fit_method = 'active-window-linear';
    end

    if use_median_speed_fallback
        if ~isfinite(fit.c) || (~isfinite(fit.R2)) || fit.R2 < R2_fallback_threshold
            fit2 = fit_speed_by_median_local_velocity(t_sol, Lfit, fit_start_eff);
            if isfinite(fit2.c)
                fit = fit2;
                fit_method = 'median-local-velocity';
            end
        end
    end

    out = struct();
    out.D = D;
    out.gamma = gamma;
    out.c = fit.c;
    out.slope = fit.slope;
    out.R2 = fit.R2;
    out.nFit = fit.nFit;
    out.fit_method = fit_method;
    out.fit_start_eff = fit_start_eff;
end

%% ========================================================================
%  REDUCED STEFAN CURVE
% ========================================================================
function stefan = compute_stefan_curve(a, ~)
    F = @(U) U .* (1 - U) .* (U - a);
    cStar = sqrt(2) * (0.5 - a);

    cVals = [linspace(1e-4, 0.98*cStar, 180), ...
        cStar - logspace(log10(0.02*cStar), -8, 120)].';
    cVals = unique(cVals, 'stable');
    cVals = cVals(cVals > 0 & cVals < cStar);

    kappaVals = nan(size(cVals));
    for i = 1:numel(cVals)
        Wfront = front_slope_from_phase_plane(F, cVals(i));
        kappaVals(i) = -cVals(i) / Wfront;
    end

    keep = isfinite(kappaVals) & kappaVals > 0 & isfinite(cVals);
    kappaVals = kappaVals(keep);
    cVals = cVals(keep);

    [kappaVals, idx] = sort(kappaVals);
    cVals = cVals(idx);

    stefan = struct();
    stefan.kappa = kappaVals;
    stefan.c = cVals;
end

function kappa = corrected_kappa(gamma, c, D)
    c2 = c.^2;
    kappa = 2 .* gamma .* c2 ./ (2 .* c2 + D);
end

function Wfront = front_slope_from_phase_plane(F, c)
    Ustart = 1 - 1e-7;
    Uend = 1e-8;
    Uspan = linspace(Ustart, Uend, 2500).';
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'MaxStep', 1e-3);

    Fp1 = (F(1) - F(1 - 1e-6)) / 1e-6;
    lambda = (-c + sqrt(c^2 - 4*Fp1)) / 2;
    Wstart = lambda * (Ustart - 1);

    [~, W] = ode15s(@(U,W) (-c.*W - F(U)) ./ W, Uspan, Wstart, opts);
    W = W(:);
    Wfront = W(end);
end

%% ========================================================================
%  FITTING
% ========================================================================
function fit = fit_interface_speed_active_window(t, L, overlap_max, overlap_tol, ...
                                                  fit_start_eff, fit_last_fraction)
    t = t(:);
    L = L(:);
    overlap_max = overlap_max(:);

    active = (t >= fit_start_eff) & ...
             isfinite(L) & ...
             isfinite(overlap_max) & ...
             (overlap_max > overlap_tol);

    idx_all = find(active);
    if numel(idx_all) < 8
        fit = empty_fit();
        return;
    end

    breaks = [0; find(diff(idx_all) > 1); numel(idx_all)];
    best_block = [];
    best_len = 0;

    for b = 1:numel(breaks)-1
        block = idx_all(breaks(b)+1:breaks(b+1));
        if numel(block) > best_len
            best_len = numel(block);
            best_block = block;
        end
    end

    if numel(best_block) < 8
        fit = empty_fit();
        return;
    end

    nBlock = numel(best_block);
    iStart = max(1, floor((1 - fit_last_fraction) * nBlock));
    fit_idx = best_block(iStart:end);

    if numel(fit_idx) < 8
        fit_idx = best_block;
    end

    fit = fit_one_window(t, L, fit_idx);
end

function fit = fit_one_window(t, L, idx)
    idx = idx(:);
    if numel(idx) < 2
        fit = empty_fit();
        return;
    end

    p = polyfit(t(idx), L(idx), 1);
    slope = p(1);
    intercept = p(2);
    c = slope;

    Lfit = polyval(p, t(idx));
    ydat = L(idx);

    SSres = sum((ydat - Lfit).^2);
    SStot = sum((ydat - mean(ydat)).^2);

    if SStot > 0
        R2 = 1 - SSres/SStot;
    else
        R2 = NaN;
    end

    fit = struct();
    fit.slope = slope;
    fit.intercept = intercept;
    fit.c = c;
    fit.R2 = R2;
    fit.nFit = numel(idx);
end

function fit = fit_speed_by_median_local_velocity(t, L, fit_start_eff)
    t = t(:);
    L = L(:);

    valid = (t >= fit_start_eff) & isfinite(L);
    tv = t(valid);
    Lv = L(valid);

    if numel(tv) < 10
        fit = empty_fit();
        return;
    end

    dL = diff(Lv);
    dt = diff(tv);
    good = isfinite(dL) & isfinite(dt) & (dt > 0);

    c_local = dL(good) ./ dt(good);
    c_local = c_local(isfinite(c_local));

    if numel(c_local) < 5
        fit = empty_fit();
        return;
    end

    n = numel(c_local);
    iStart = max(1, floor(0.5*n));
    c_late = c_local(iStart:end);
    c_med = median(c_late, 'omitnan');

    fit = struct();
    fit.slope = c_med;
    fit.intercept = NaN;
    fit.c = c_med;
    fit.R2 = NaN;
    fit.nFit = numel(c_late);
end

function fit = empty_fit()
    fit = struct();
    fit.slope = NaN;
    fit.intercept = NaN;
    fit.c = NaN;
    fit.R2 = NaN;
    fit.nFit = 0;
end

%% ========================================================================
%  PDE HELPERS
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

function Jpat = jacobian_pattern(Nx, Lmat)
    P = spones(Lmat);
    I = speye(Nx);
    Jpat = [P + I, I;
            I,     P + I];
end

function L = neumann_laplacian_1d(Nx, dx)
    e = ones(Nx,1);
    L = spdiags([e -2*e e], -1:1, Nx, Nx);
    L(1,2) = 2;
    L(Nx,Nx-1) = 2;
    L = L / dx^2;
end

function [u0, v0] = initial_condition_gap(x, A, w, x_uedge, x_vedge)
    H1 = 0.5 * (1 + tanh((x - x_uedge)/w));
    H2 = 0.5 * (1 + tanh((x - x_vedge)/w));
    u0 = A * (1 - H1);
    v0 = A * H2;
end

function [u0, v0] = initial_condition_stefan_small_branch(x, A, w, L0, a1, D)
    U = linspace(0, 1-1e-8, 5000);
    R = U .* (1-U) .* (U-a1);
    I = cumtrapz(U, R);
    totalI = I(end);
    tailI = max(totalI - I, 0);
    W0 = -sqrt(2*tailI);

    invW = 1 ./ W0;
    invW(~isfinite(invW)) = 0;
    z = cumtrapz(U, invW);

    z = z(:);
    U = U(:);
    [zUnique, ia] = unique(z, 'stable');
    UUnique = U(ia);

    zz = x(:) - L0;
    u0 = interp1(zUnique, UUnique, zz, 'linear', 1);
    u0(zz >= 0) = 0;
    u0 = max(0, min(A, A*u0));

    vWidth = max(w, 4*sqrt(max(D, eps)));
    H = 0.5 * (1 + tanh(zz / vWidth));
    v0 = A * H;
    v0 = max(0, min(A, v0));
end

function xpeak = subgrid_peak_location(x, w, idxMax)
    xpeak = x(idxMax);

    if idxMax <= 1 || idxMax >= numel(w)
        return;
    end

    wLocal = w([idxMax-1, idxMax, idxMax+1]);
    if any(~isfinite(wLocal)) || any(wLocal <= 0)
        return;
    end

    g = log(wLocal);
    denom = g(1) - 2*g(2) + g(3);
    if abs(denom) < eps
        return;
    end

    offset = 0.5 * (g(1) - g(3)) / denom;
    offset = max(min(offset, 1), -1);
    dx = x(2) - x(1);
    xpeak = x(idxMax) + offset * dx;
end

function val = pick(s, name, defaultVal)
    if isfield(s, name)
        val = s.(name);
    else
        val = defaultVal;
    end
end
