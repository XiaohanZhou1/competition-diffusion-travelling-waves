function Ode15sJacobian()
% Ode15sJacobian_surface_contour
%
% Solve
%   u_t = u_xx + u(1-u)(u-a1) - uv/delta
%   v_t = D v_xx + v(1-v) - gamma*uv/delta
%
% on [0, Xmax] with homogeneous Neumann BCs,
% using method of lines + ode15s + sparse Jacobian.
%
% Then compute signed wave speed
%   c = - dL/dt
% where L(t) is by default taken as
%   L(t) = argmax_x (u(x,t)v(x,t)).
%
% Only two figures are produced:
%   1) contour plot of c(D,gamma)
%   2) surface plot of c(D,gamma)

    clearvars -except ans;
    clc;
    close all;

    %% ============================================================
    %  USER PARAMETERS
    % ============================================================

    % -------- parameter scan --------
    D_list     = linspace(0.2, 10, 16);
    gamma_list = linspace(0.05, 5, 21);

    % -------- numerical options --------
    opts = struct();

    opts.delta = 1e-3;
    opts.a1    = 0.1;
    opts.A     = 1.0;

    opts.Xmax = 250;
    opts.Nx   = 500;      % if too slow, reduce to 350 first

    opts.tEnd  = 180;     % slightly longer than before
    opts.NtOut = 280;     % number of output times

    opts.gap = 80;
    opts.wIC = 1.0;

    % threshold for identifying fronts (diagnostic only)
    opts.eta_u = 0.01;
    opts.eta_v = 0.01;

    % tolerance used in u=v interface finder
    opts.active_tol_equal = 1e-10;

    % overlap threshold for saying the two species really interact
    opts.overlap_tol = 1e-8;

    % fit settings
    opts.t_fit_start = 10;
    opts.fit_last_fraction = 0.40;

    % choose interface for wave speed:
    %   'uvmax'  (recommended, robust)
    %   'equal'  (use u=v)
    opts.interface_for_speed = 'uvmax';

    % ode15s options
    opts.RelTol  = 1e-5;
    opts.AbsTol  = 1e-7;
    opts.MaxStep = 1.0;

    opts.verbose = false;

    %% ============================================================
    %  STORAGE
    % ============================================================

    ND     = numel(D_list);
    Ngamma = numel(gamma_list);

    c_mat     = nan(Ngamma, ND);
    slope_mat = nan(Ngamma, ND);
    R2_mat    = nan(Ngamma, ND);

    fprintf('Starting parameter sweep using ode15s + sparse Jacobian...\n');
    fprintf('Number of D values     = %d\n', ND);
    fprintf('Number of gamma values = %d\n', Ngamma);
    fprintf('Total ODE solves       = %d\n\n', ND * Ngamma);

    tStartAll = tic;

    for j = 1:ND
        D = D_list(j);

        for i = 1:Ngamma
            gamma = gamma_list(i);
            case_id = (j-1)*Ngamma + i;

            fprintf('Case (%d/%d): D = %.4f, gamma = %.4f\n', ...
                case_id, ND*Ngamma, D, gamma);

            try
                out = run_one_case_ode15s(D, gamma, opts);

                c_mat(i,j)     = out.c;
                slope_mat(i,j) = out.slope;
                R2_mat(i,j)    = out.R2;

                fprintf('    c = %+ .6e,  R2 = %.4f,  fit points = %d\n', ...
                    out.c, out.R2, out.nFit);

            catch ME
                warning('Failed at D = %.6f, gamma = %.6f: %s', ...
                    D, gamma, ME.message);
            end
        end
    end

    elapsed = toc(tStartAll);

    fprintf('\nSweep finished.\n');
    fprintf('Elapsed time = %.2f seconds.\n', elapsed);

    %% ============================================================
    %  BASIC SUMMARY
    % ============================================================

    nFinite = nnz(isfinite(c_mat));
    fprintf('Finite c values found: %d / %d\n', nFinite, numel(c_mat));

    if nFinite == 0
        error(['No finite wave speeds were obtained. ', ...
               'Try reducing gap, increasing tEnd, or lowering overlap_tol.']);
    end

    %% ============================================================
    %  MESHGRID
    % ============================================================

    [DD, GG] = meshgrid(D_list, gamma_list);

    %% ============================================================
    %  CONTOUR PLOT
    % ============================================================

    c_plot = fill_nan_for_plot(c_mat);

    finite_vals = c_plot(isfinite(c_plot));
    cmin = min(finite_vals);
    cmax = max(finite_vals);

    dc_fill = 0.1;
    low_level  = dc_fill * floor(cmin / dc_fill);
    high_level = dc_fill * ceil(cmax / dc_fill);

    if abs(high_level - low_level) < 1e-12
        levels_fill = linspace(cmin - 0.5*dc_fill, cmax + 0.5*dc_fill, 12);
    else
        levels_fill = low_level:dc_fill:high_level;
    end

    levels_fill = levels_fill(isfinite(levels_fill));
    levels_fill = unique(levels_fill);

    if numel(levels_fill) < 2
        levels_fill = linspace(cmin, cmax, 12);
    end

    figure;
    contourf(DD, GG, c_plot, levels_fill, 'LineColor', 'none');
    colorbar;
    hold on;

    if cmin <= 0 && cmax >= 0
        contour(DD, GG, c_plot, [0 0], 'k', 'LineWidth', 2.0);
    end

    xlabel('$D$', 'Interpreter', 'latex');
    ylabel('$\gamma$', 'Interpreter', 'latex');
    box on;
    grid off;

    %% ============================================================
    %  SURFACE PLOT
    % ============================================================

    figure;
    surf(DD, GG, c_mat, 'EdgeColor', 'none');
    colorbar;
    xlabel('$D$', 'Interpreter', 'latex');
    ylabel('$\gamma$', 'Interpreter', 'latex');
    zlabel('$c$', 'Interpreter', 'latex');
    view(135, 30);
    zlim([-0.5 3]);
    box on;

    fprintf('\nDone.\n');

end

%% ========================================================================
%  SOLVE ONE PARAMETER CASE
% ========================================================================
function out = run_one_case_ode15s(D, gamma, opts)

    % ----- parameters -----
    delta = pick(opts, 'delta', 1e-3);
    a1    = pick(opts, 'a1', 0.1);
    A     = pick(opts, 'A', 1.0);

    Xmax = pick(opts, 'Xmax', 250);
    Nx   = pick(opts, 'Nx', 500);

    tEnd  = pick(opts, 'tEnd', 180);
    NtOut = pick(opts, 'NtOut', 280);

    gap = pick(opts, 'gap', 80);
    wIC = pick(opts, 'wIC', 1.0);

    active_tol_equal = pick(opts, 'active_tol_equal', 1e-10);
    overlap_tol      = pick(opts, 'overlap_tol', 1e-8);

    t_fit_start       = pick(opts, 't_fit_start', 10);
    fit_last_fraction = pick(opts, 'fit_last_fraction', 0.40);

    interface_for_speed = pick(opts, 'interface_for_speed', 'uvmax');

    RelTol  = pick(opts, 'RelTol', 1e-5);
    AbsTol  = pick(opts, 'AbsTol', 1e-7);
    MaxStep = pick(opts, 'MaxStep', 1.0);

    verbose = pick(opts, 'verbose', false);

    % ----- grid -----
    x  = linspace(0, Xmax, Nx).';
    dx = x(2) - x(1);

    Lmat = neumann_laplacian_1d(Nx, dx);

    % ----- initial condition -----
    xmid    = 0.5 * Xmax;
    x_uedge = xmid - gap/2;
    x_vedge = xmid + gap/2;

    [u0, v0] = initial_condition_gap(x, A, wIC, x_uedge, x_vedge);
    y0 = [u0; v0];

    % ----- output times -----
    tspan = linspace(0, tEnd, NtOut);

    % ----- ODE options -----
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

    % ----- solve -----
    [t_sol, y_sol] = ode15s(rhs, tspan, y0, ode_opts);

    Nt = numel(t_sol);

    u = y_sol(:, 1:Nx);
    v = y_sol(:, Nx+1:2*Nx);

    % ----- compute interfaces for ALL times -----
    overlap_max = nan(Nt,1);
    L_uvmax     = nan(Nt,1);
    L_equal     = nan(Nt,1);

    prevL = NaN;

    for k = 1:Nt
        uk = u(k,:).';
        vk = v(k,:).';

        % overlap-based interface
        wk = uk .* vk;
        [wmax, idxMax] = max(wk);
        overlap_max(k) = wmax;

        if wmax > overlap_tol
            L_uvmax(k) = subgrid_peak_location(x, wk, idxMax);
        end

        % u = v interface
        Lk = find_interface_uv_equal(x, uk, vk, prevL, active_tol_equal);

        if ~isnan(Lk)
            L_equal(k) = Lk;
            prevL = Lk;
        end
    end

    % ----- choose fit start based on interaction -----
    idxInteract = find(overlap_max > overlap_tol, 1, 'first');

    if isempty(idxInteract)
        fit_start_eff = t_fit_start;
    else
        fit_start_eff = max(t_fit_start, t_sol(idxInteract));
    end

    % ----- choose interface for speed -----
    switch lower(interface_for_speed)
        case 'uvmax'
            Lfit = L_uvmax;
        case 'equal'
            Lfit = L_equal;
        otherwise
            error('Unknown interface_for_speed = %s', interface_for_speed);
    end

    % ----- fit speed -----
    fit = fit_interface_speed(t_sol, Lfit, overlap_max, overlap_tol, ...
                              fit_start_eff, fit_last_fraction);

    if verbose
        fprintf('D = %.4f, gamma = %.4f, c = %.6e, R2 = %.6f\n', ...
            D, gamma, fit.c, fit.R2);
    end

    out = struct();
    out.D = D;
    out.gamma = gamma;
    out.c = fit.c;
    out.slope = fit.slope;
    out.R2 = fit.R2;
    out.nFit = fit.nFit;

    out.t = t_sol;
    out.x = x;
    out.L_uvmax = L_uvmax;
    out.L_equal = L_equal;
    out.overlap_max = overlap_max;
    out.fit_start_eff = fit_start_eff;
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
%  JACOBIAN
% ========================================================================
function J = jacobian_fun(y, Lmat, D, delta, gamma, a1)

    Nx = numel(y)/2;

    u = y(1:Nx);
    v = y(Nx+1:end);

    % R(u)=u(1-u)(u-a1) = -u^3 + (1+a1)u^2 - a1*u
    Rp = -3*u.^2 + 2*(1+a1)*u - a1;

    % S(v)=v(1-v)
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

    L(1,2)     = 2;
    L(Nx,Nx-1) = 2;

    L = L / dx^2;
end

%% ========================================================================
%  INITIAL CONDITION
% ========================================================================
function [u0, v0] = initial_condition_gap(x, A, w, x_uedge, x_vedge)

    H1 = 0.5 * (1 + tanh((x - x_uedge)/w));
    H2 = 0.5 * (1 + tanh((x - x_vedge)/w));

    u0 = A * (1 - H1);
    v0 = A * H2;
end

%% ========================================================================
%  FIT INTERFACE SPEED
% ========================================================================
function fit = fit_interface_speed(t, L, overlap_max, overlap_tol, ...
                                   fit_start_eff, fit_last_fraction)

    t = t(:);
    L = L(:);
    overlap_max = overlap_max(:);

    valid_idx = find( (t >= fit_start_eff) & ...
                      isfinite(L) & ...
                      isfinite(overlap_max) & ...
                      (overlap_max > overlap_tol) );

    if numel(valid_idx) >= 8
        nValid = numel(valid_idx);
        iStart = max(1, floor((1 - fit_last_fraction)*nValid));
        fit_idx = valid_idx(iStart:end);

        p = polyfit(t(fit_idx), L(fit_idx), 1);

        slope = p(1);
        intercept = p(2);
        c = -slope;

        Lfit = polyval(p, t(fit_idx));
        ydat = L(fit_idx);

        SSres = sum((ydat - Lfit).^2);
        SStot = sum((ydat - mean(ydat)).^2);

        if SStot > 0
            R2 = 1 - SSres/SStot;
        else
            R2 = NaN;
        end

        nFit = numel(fit_idx);

    elseif numel(valid_idx) >= 2
        fit_idx = valid_idx;

        p = polyfit(t(fit_idx), L(fit_idx), 1);

        slope = p(1);
        intercept = p(2);
        c = -slope;

        Lfit = polyval(p, t(fit_idx));
        ydat = L(fit_idx);

        SSres = sum((ydat - Lfit).^2);
        SStot = sum((ydat - mean(ydat)).^2);

        if SStot > 0
            R2 = 1 - SSres/SStot;
        else
            R2 = NaN;
        end

        nFit = numel(fit_idx);

    else
        slope = NaN;
        intercept = NaN;
        c = NaN;
        R2 = NaN;
        nFit = 0;
    end

    fit = struct();
    fit.slope = slope;
    fit.intercept = intercept;
    fit.c = c;
    fit.R2 = R2;
    fit.nFit = nFit;
end

%% ========================================================================
%  u = v INTERFACE
% ========================================================================
function Lx = find_interface_uv_equal(x, urow, vrow, prevL, active_tol)

    d = urow - vrow;
    idx = find(d(1:end-1).*d(2:end) <= 0);

    xcand = [];
    wcand = [];

    for j = 1:numel(idx)
        i = idx(j);

        d1 = d(i);
        d2 = d(i+1);

        if abs(d2 - d1) < 1e-14
            theta = 0.5;
        else
            theta = -d1 / (d2 - d1);
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

    if ~isempty(xcand)
        if ~isnan(prevL)
            [~, ibest] = min(abs(xcand - prevL));
        else
            [~, ibest] = max(wcand);
        end
        Lx = xcand(ibest);
        return;
    end

    % fallback: minimal |u-v|
    absd = abs(d);

    if ~isnan(prevL)
        window = 20;
        [~, ic] = min(abs(x - prevL));
        i1 = max(1, ic-window);
        i2 = min(numel(x), ic+window);

        [~, iloc] = min(absd(i1:i2));
        iBest = i1 + iloc - 1;
    else
        [~, iBest] = min(absd);
    end

    ui = urow(iBest);
    vi = vrow(iBest);

    if min(ui,vi) < active_tol
        Lx = NaN;
    else
        Lx = x(iBest);
    end
end

%% ========================================================================
%  SUBGRID LOCATION OF MAXIMUM OF w = u*v
% ========================================================================
function xpeak = subgrid_peak_location(x, w, idxMax)

    n = numel(x);

    if idxMax <= 1 || idxMax >= n
        xpeak = x(idxMax);
        return;
    end

    xx = x(idxMax-1:idxMax+1);
    ww = w(idxMax-1:idxMax+1);

    % quadratic fit
    p = polyfit(xx, ww, 2);

    if abs(p(1)) < 1e-14
        xpeak = x(idxMax);
        return;
    end

    xv = -p(2)/(2*p(1));

    % keep vertex inside the local interval
    xmin = xx(1);
    xmax = xx(3);

    if xv < xmin || xv > xmax || ~isfinite(xv)
        xpeak = x(idxMax);
    else
        xpeak = xv;
    end
end

%% ========================================================================
%  FILL NaNs FOR PLOTTING
% ========================================================================
function Afill = fill_nan_for_plot(A)

    Afill = A;

    finite_mask = isfinite(Afill);

    if nnz(finite_mask) == 0
        error('fill_nan_for_plot: the matrix contains no finite values.');
    end

    if nnz(finite_mask) == 1
        only_val = Afill(finite_mask);
        Afill(~finite_mask) = only_val;
        return;
    end

    Afill = fillmissing(Afill, 'linear', 1, 'EndValues', 'nearest');
    Afill = fillmissing(Afill, 'linear', 2, 'EndValues', 'nearest');

    finite_mask = isfinite(Afill);
    if any(~finite_mask(:))
        m = mean(Afill(finite_mask), 'omitnan');
        Afill(~finite_mask) = m;
    end
end

%% ========================================================================
%  PICK FIELD OR DEFAULT
% ========================================================================
function val = pick(s, name, defaultVal)

    if isfield(s, name)
        val = s.(name);
    else
        val = defaultVal;
    end
end