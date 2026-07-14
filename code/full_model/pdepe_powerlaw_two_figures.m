%% ============================================================
%  pdepe_powerlaw_two_figures.m
%
%  Use PDEPE to solve the two-species competition-diffusion system:
%
%    u_t = u_xx + u(1-u)(u-a1) - uv/delta,
%    v_t = D v_xx + v(1-v)     - gamma uv/delta.
%
%  The interface L(t) is defined by u = v.
%  The signed wave speed is defined by
%
%       c = - dL/dt.
%
%  Then extract gamma_crit(D) from c = 0 and plot:
%    Figure 1: gamma_crit vs D for D < 1;
%    Figure 2: gamma_crit vs D for 1 < D < 10.
%% ============================================================

clear;
clc;
close all;

%% ============================================================
%  PARAMETER GRID
%% ============================================================

% D values for small-D and large-D regimes
D_list = [0.1:0.1:1.0, 1.2:0.4:9.2];

% gamma values. Increase Ngamma for smoother/more accurate critical curve.
gamma_list = linspace(0.1, 5.2, 38);

ND     = numel(D_list);
Ngamma = numel(gamma_list);

%% ============================================================
%  NUMERICAL OPTIONS
%% ============================================================

opts = struct();

opts.delta = 1e-3;
opts.a1    = 0.1;
opts.A     = 1.0;

opts.Xmax = 250;
opts.Nx   = 600;

opts.tEnd = 100;
opts.Nt   = 280;

opts.gap = 80;
opts.wIC = 1.0;

opts.eta_u = 0.01;
opts.eta_v = 0.01;

opts.active_tol  = 1e-5;
opts.t_fit_start = 20;

opts.verbose = false;

%% ============================================================
%  STORAGE
%% ============================================================

c_mat     = nan(Ngamma, ND);
slope_mat = nan(Ngamma, ND);
R2_mat    = nan(Ngamma, ND);
tc_mat    = nan(Ngamma, ND);

fprintf('Starting PDEPE parameter sweep...\n');
fprintf('Number of D values     = %d\n', ND);
fprintf('Number of gamma values = %d\n', Ngamma);
fprintf('Total PDE solves       = %d\n\n', ND * Ngamma);

tic;

%% ============================================================
%  PARAMETER SWEEP
%% ============================================================

for j = 1:ND
    D = D_list(j);

    for i = 1:Ngamma
        gamma = gamma_list(i);

        case_id = (j-1)*Ngamma + i;

        fprintf('Case %d/%d: D = %.4f, gamma = %.4f\n', ...
            case_id, ND*Ngamma, D, gamma);

        try
            [c_est, result] = compute_c_from_pdepe_signed(D, gamma, opts);

            c_mat(i,j)     = c_est;
            slope_mat(i,j) = result.slopeL;
            R2_mat(i,j)    = result.R2;
            tc_mat(i,j)    = result.tc;

            fprintf('    c = %+ .6e, R2 = %.4f, tc = %.4f\n', ...
                c_est, result.R2, result.tc);

        catch ME
            warning('Failed at D = %.4f, gamma = %.4f: %s', ...
                D, gamma, ME.message);
        end
    end
end

fprintf('\nSweep finished in %.2f seconds.\n', toc);
fprintf('Finite c values: %d / %d\n', nnz(isfinite(c_mat)), numel(c_mat));

%% ============================================================
%  ESTIMATE CRITICAL CURVE gamma_crit(D)
%% ============================================================

gamma_crit = nan(1, ND);

for j = 1:ND
    gamma_crit(j) = estimate_gamma_critical(gamma_list, c_mat(:,j));
end

mask_valid = isfinite(gamma_crit) & gamma_crit > 0 & D_list > 0;

fprintf('\nValid critical points: %d / %d\n', nnz(mask_valid), ND);

%% ============================================================
%  FIGURE 1: SMALL-D POWER-LAW FIT, D < 1
%% ============================================================

blueColor = [0.0000, 0.4470, 0.7410];

mask_small = mask_valid & (D_list < 1);

D_small = D_list(mask_small);
g_small = gamma_crit(mask_small);

figure('Color','w');
hold on;

if numel(D_small) >= 2

    p_small = polyfit(log(D_small(:)), log(g_small(:)), 1);

    alpha_small = p_small(1);
    K_small     = exp(p_small(2));

    D_small_plot = logspace(log10(min(D_small)), log10(max(D_small)), 300);
    g_small_fit  = K_small * D_small_plot.^alpha_small;

    loglog(D_small, g_small, 'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', 6, ...
        'DisplayName', '$\gamma(D)$ data');

    loglog(D_small_plot, g_small_fit, '-', ...
        'Color', blueColor, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('$%.4fD^{%.4f}$', K_small, alpha_small));

    fprintf('\n===== SMALL-D POWER-LAW FIT: D < 1 =====\n');
    fprintf('gamma_crit(D) = %.8f D^{%.8f}\n', K_small, alpha_small);

else
    warning('Not enough valid points for small-D fit.');
end

xlabel('$D$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\gamma$', 'Interpreter', 'latex', 'FontSize', 20);

legend('Location', 'best', ...
    'Interpreter', 'latex', ...
    'FontSize', 11);

set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 16, ...
    'LineWidth', 1.1);

grid on;
box on;

%% ============================================================
%  FIGURE 2: LARGE-D POWER-LAW FIT, 1 < D < 10
%% ============================================================

mask_large = mask_valid & (D_list > 1) & (D_list < 10);

D_large = D_list(mask_large);
g_large = gamma_crit(mask_large);

figure('Color','w');
hold on;

if numel(D_large) >= 2

    p_large = polyfit(log(D_large(:)), log(g_large(:)), 1);

    alpha_large = p_large(1);
    K_large     = exp(p_large(2));

    D_large_plot = logspace(log10(min(D_large)), log10(max(D_large)), 300);
    g_large_fit  = K_large * D_large_plot.^alpha_large;

    loglog(D_large, g_large, 'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', 6, ...
        'DisplayName', '$\gamma(D)$ data');

    loglog(D_large_plot, g_large_fit, '-', ...
        'Color', blueColor, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('$%.4fD^{%.4f}$', K_large, alpha_large));

    fprintf('\n===== LARGE-D POWER-LAW FIT: 1 < D < 10 =====\n');
    fprintf('gamma_crit(D) = %.8f D^{%.8f}\n', K_large, alpha_large);

else
    warning('Not enough valid points for large-D fit.');
end

xlabel('$D$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\gamma$', 'Interpreter', 'latex', 'FontSize', 20);

legend('Location', 'best', ...
    'Interpreter', 'latex', ...
    'FontSize', 11);

set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 16, ...
    'LineWidth', 1.1);

grid on;
box on;

%% ============================================================
%  LOCAL FUNCTIONS
%% ============================================================

function [c_est, result] = compute_c_from_pdepe_signed(D, gamma, opts)
% Solve the coupled competition PDE using PDEPE and estimate signed speed.
%
% Interface:
%   L(t) is defined by u = v.
%
% Speed convention:
%   L(t) ~ alpha + slope * t,
%   c = -slope.
%
% Hence:
%   c > 0 means L(t) moves left;
%   c < 0 means L(t) moves right.
% Adjust this interpretation if your paper uses the opposite convention.

    delta       = pick(opts, 'delta',       1e-3);
    a1          = pick(opts, 'a1',          0.1);
    A           = pick(opts, 'A',           1.0);

    Xmax        = pick(opts, 'Xmax',        250);
    Nx          = pick(opts, 'Nx',          600);

    tEnd        = pick(opts, 'tEnd',        80);
    Nt          = pick(opts, 'Nt',          250);

    gap         = pick(opts, 'gap',         80);
    wIC         = pick(opts, 'wIC',         1.0);

    eta_u       = pick(opts, 'eta_u',       0.01);
    eta_v       = pick(opts, 'eta_v',       0.01);

    active_tol  = pick(opts, 'active_tol',  1e-5);
    t_fit_start = pick(opts, 't_fit_start', 20);

    verbose     = pick(opts, 'verbose',     false);

    R = @(u) u.*(1-u).*(u-a1);
    S = @(v) v.*(1-v);

    x    = linspace(0, Xmax, Nx);
    t    = linspace(0, tEnd, Nt);
    tcol = t(:);
    m    = 0;

    xmid    = 0.5 * Xmax;
    x_uedge = xmid - gap/2;
    x_vedge = xmid + gap/2;

    sol = pdepe(m, ...
        @(x,t,U,DUdx) pdefun(x,t,U,DUdx,D,delta,gamma,R,S), ...
        @(x) icfun_gap(x,A,wIC,x_uedge,x_vedge), ...
        @bcfun_neumann, ...
        x, t);

    u = sol(:,:,1);
    v = sol(:,:,2);

    Xu = nan(Nt,1);
    Xv = nan(Nt,1);

    for k = 1:Nt
        idx_u = find(u(k,:) >= eta_u, 1, 'last');
        if ~isempty(idx_u)
            Xu(k) = x(idx_u);
        end

        idx_v = find(v(k,:) >= eta_v, 1, 'first');
        if ~isempty(idx_v)
            Xv(k) = x(idx_v);
        end
    end

    tc_idx = find(Xu >= Xv, 1, 'first');

    if isempty(tc_idx)
        tc = NaN;
    else
        tc = t(tc_idx);
    end

    L = nan(Nt,1);

    if ~isempty(tc_idx)
        prevL = NaN;

        for k = tc_idx:Nt
            Lk = find_interface_uv_equal(x, u(k,:), v(k,:), prevL, active_tol);
            L(k) = Lk;

            if ~isnan(Lk)
                prevL = Lk;
            end
        end
    end

    if isnan(tc)
        fit_start_eff = t_fit_start;
    else
        fit_start_eff = max(t_fit_start, tc);
    end

    fit_mask = (tcol > fit_start_eff) & isfinite(L);

    if nnz(fit_mask) >= 8

        % Fit only the last 40% of valid points.
        valid_idx = find(fit_mask);
        nValid = numel(valid_idx);
        iStart = max(1, floor(0.60*nValid));
        fit_idx = valid_idx(iStart:end);

        p = polyfit(tcol(fit_idx), L(fit_idx), 1);

        slopeL     = p(1);
        interceptL = p(2);

        c_est = -slopeL;

        L_fit = polyval(p, tcol(fit_idx));
        ydata = L(fit_idx);

        SSres = sum((ydata - L_fit).^2);
        SStot = sum((ydata - mean(ydata)).^2);

        if SStot > 0
            R2 = 1 - SSres/SStot;
        else
            R2 = NaN;
        end

    elseif nnz(fit_mask) >= 2

        fit_idx = find(fit_mask);

        p = polyfit(tcol(fit_idx), L(fit_idx), 1);

        slopeL     = p(1);
        interceptL = p(2);

        c_est = -slopeL;

        L_fit = polyval(p, tcol(fit_idx));
        ydata = L(fit_idx);

        SSres = sum((ydata - L_fit).^2);
        SStot = sum((ydata - mean(ydata)).^2);

        if SStot > 0
            R2 = 1 - SSres/SStot;
        else
            R2 = NaN;
        end

    else
        slopeL     = NaN;
        interceptL = NaN;
        c_est      = NaN;
        R2         = NaN;
    end

    if verbose
        fprintf('D = %.6g, gamma = %.6g, c = %.8f, tc = %.8f, R2 = %.6f\n', ...
            D, gamma, c_est, tc, R2);
    end

    result = struct();

    result.D          = D;
    result.gamma      = gamma;
    result.c          = c_est;
    result.tc         = tc;
    result.x          = x;
    result.t          = tcol;
    result.u          = u;
    result.v          = v;
    result.Xu         = Xu;
    result.Xv         = Xv;
    result.L          = L;
    result.slopeL     = slopeL;
    result.interceptL = interceptL;
    result.R2         = R2;
    result.fit_mask   = fit_mask;
    result.fit_start  = fit_start_eff;
end

function val = pick(s, name, defaultVal)

    if isfield(s, name)
        val = s.(name);
    else
        val = defaultVal;
    end
end

function [c,f,s] = pdefun(~,~,U,DUdx,D,delta,gamma,R,S)

    u  = U(1);
    v  = U(2);

    ux = DUdx(1);
    vx = DUdx(2);

    c = [1; 1];

    f = [ux; D*vx];

    s = [R(u) - (u*v)/delta;
         S(v) - gamma*(u*v)/delta];
end

function U0 = icfun_gap(x,A,w,x_uedge,x_vedge)

    H1 = 0.5*(1 + tanh((x - x_uedge)/w));
    H2 = 0.5*(1 + tanh((x - x_vedge)/w));

    u0 = A*(1 - H1);
    v0 = A*H2;

    U0 = [u0; v0];
end

function [pl,ql,pr,qr] = bcfun_neumann(~,~,~,~,~)

    pl = [0; 0];
    ql = [1; 1];

    pr = [0; 0];
    qr = [1; 1];
end

function Lx = find_interface_uv_equal(x, urow, vrow, prevL, active_tol)

    d = urow - vrow;

    idx = find(d(1:end-1).*d(2:end) <= 0);

    if isempty(idx)
        Lx = NaN;
        return;
    end

    xcand = [];
    wcand = [];

    for j = 1:numel(idx)
        i = idx(j);

        d1 = d(i);
        d2 = d(i+1);

        if abs(d2 - d1) < 1e-14
            theta = 0.5;
        else
            theta = -d1/(d2 - d1);
        end

        if theta < 0 || theta > 1
            continue;
        end

        xi = x(i) + theta*(x(i+1)-x(i));

        ui = urow(i) + theta*(urow(i+1)-urow(i));
        vi = vrow(i) + theta*(vrow(i+1)-vrow(i));

        if min(ui,vi) < active_tol
            continue;
        end

        xcand(end+1,1) = xi; %#ok<AGROW>
        wcand(end+1,1) = ui + vi; %#ok<AGROW>
    end

    if isempty(xcand)
        Lx = NaN;
        return;
    end

    if ~isnan(prevL)
        [~, jbest] = min(abs(xcand - prevL));
    else
        [~, jbest] = max(wcand);
    end

    Lx = xcand(jbest);
end

function gamma_crit = estimate_gamma_critical(gamma_list, c_col)

    gamma_crit = NaN;

    valid = isfinite(c_col);

    if nnz(valid) < 2
        return;
    end

    g = gamma_list(valid);
    c = c_col(valid);

    [cmin_abs, idxmin] = min(abs(c));

    if cmin_abs < 1e-12
        gamma_crit = g(idxmin);
        return;
    end

    idx = find(c(1:end-1).*c(2:end) < 0);

    if isempty(idx)
        return;
    end

    score = inf(size(idx));

    for k = 1:numel(idx)
        ii = idx(k);
        score(k) = min(abs([c(ii), c(ii+1)]));
    end

    [~, bestk] = min(score);
    ii = idx(bestk);

    g1 = g(ii);
    g2 = g(ii+1);
    c1 = c(ii);
    c2 = c(ii+1);

    gamma_crit = g1 - c1*(g2 - g1)/(c2 - c1);
end