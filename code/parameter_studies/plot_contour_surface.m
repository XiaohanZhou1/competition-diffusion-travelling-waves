function plot_contour_surface()
% Solve
%   u_t = u_xx + u(1-u)(u-a1) - uv/delta
%   v_t = D v_xx + v(1-v) - gamma*uv/delta
%
% on [0, Xmax] with homogeneous Neumann boundary conditions,
% using method of lines + ode15s + sparse Jacobian.
%
% Interface options:
%   'uvmax'         : L(t) = argmax_x u v
%   'localcentroid' : local uv-weighted centroid near argmax_x uv
%   'centroid'      : global uv-weighted centroid
%   'equal'         : u = v
%
% Signed wave speed:
%   c = - dL/dt.
%
% Special treatment:
%   For D > 5 and gamma < 0.01, the simulation uses a shorter time window
%   to extract the speed during the active-overlap phase before max_x(uv)
%   collapses to numerical noise.

    clearvars -except ans;
    clc;
    close all;

    %% ============================================================
    %  USER PARAMETERS
    % ============================================================

    % -------- parameter scan --------
    D_list = linspace(0.1, 10, 31);

    % Include gamma < 0.01 so the special large-D/small-gamma branch is used.
    gamma_list = unique([ ...
        linspace(0.005, 0.10, 8), ...
        linspace(0.15, 1.00, 10), ...
        linspace(1.25, 5.00, 16) ...
    ]);

    % -------- numerical options --------
    opts = struct();

    opts.delta = 1e-3;
    opts.a1    = 0.3;
    opts.A     = 1.0;

    opts.Xmax = 250;
    opts.Nx   = 500;

    opts.tEnd  = 180;
    opts.NtOut = 280;

    opts.gap = 80;
    opts.wIC = 1.0;

    opts.eta_u = 0.01;
    opts.eta_v = 0.01;

    opts.active_tol_equal = 1e-10;
    opts.overlap_tol      = 1e-8;

    opts.t_fit_start       = 10;
    opts.fit_last_fraction = 0.40;

    opts.interface_for_speed = 'localcentroid';
    opts.local_centroid_half_window = 12;

    opts.RelTol  = 1e-5;
    opts.AbsTol  = 1e-7;
    opts.MaxStep = 1.0;

    opts.verbose = false;

    opts.use_median_speed_fallback = true;
    opts.R2_fallback_threshold = 0.95;

    % Special settings for very difficult region
    opts.use_early_window_for_largeD_smallgamma = true;
    opts.early_D_min = 5.0;
    opts.early_gamma_max = 0.01;

    %% ============================================================
    %  STORAGE
    % ============================================================

    ND     = numel(D_list);
    Ngamma = numel(gamma_list);

    c_mat     = nan(Ngamma, ND);
    slope_mat = nan(Ngamma, ND);
    R2_mat    = nan(Ngamma, ND);
    nFit_mat  = nan(Ngamma, ND);

    fprintf('Starting parameter sweep using ode15s + sparse Jacobian...\n');
    fprintf('Interface definition   = %s\n', opts.interface_for_speed);
    fprintf('Number of D values     = %d\n', ND);
    fprintf('Number of gamma values = %d\n', Ngamma);
    fprintf('Total ODE solves       = %d\n\n', ND * Ngamma);

    tStartAll = tic;

    %% ============================================================
    %  PARAMETER SWEEP
    % ============================================================

    for j = 1:ND
        D = D_list(j);

        for i = 1:Ngamma
            gamma = gamma_list(i);
            case_id = (j-1)*Ngamma + i;

            fprintf('Case (%d/%d): D = %.4f, gamma = %.4f\n', ...
                case_id, ND*Ngamma, D, gamma);

            try
                opts_case = opts;

                % ------------------------------------------------------------
                % Special early-window treatment for D > 5, gamma < 0.01.
                % Based on diagnostics, max_x(uv) is only significant during
                % an early active-overlap window. Running to very large t
                % mainly adds numerical tail/noise and does not improve c.
                % ------------------------------------------------------------
                if opts.use_early_window_for_largeD_smallgamma && ...
                        (D > opts.early_D_min) && ...
                        (gamma < opts.early_gamma_max)

                    fprintf('    Using early active-overlap window settings...\n');

                    opts_case.Xmax  = 400;
                    opts_case.Nx    = 800;

                    opts_case.tEnd  = 30;
                    opts_case.NtOut = 300;

                    opts_case.gap = 50;

                    opts_case.overlap_tol = 1e-13;
                    opts_case.MaxStep     = 0.2;
                    opts_case.t_fit_start = 10;

                    opts_case.local_centroid_half_window = 16;
                    opts_case.interface_for_speed = 'localcentroid';
                end

                out = run_one_case_ode15s(D, gamma, opts_case);

                c_mat(i,j)     = out.c;
                slope_mat(i,j) = out.slope;
                R2_mat(i,j)    = out.R2;
                nFit_mat(i,j)  = out.nFit;

                fprintf('    c = %+ .6e,  R2 = %.4f,  fit points = %d,  method = %s\n', ...
                    out.c, out.R2, out.nFit, out.fit_method);

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
    %  SECOND PASS: RERUN MISSING POINTS WITH AGGRESSIVE SETTINGS
    % ============================================================

    missing_mask = ~isfinite(c_mat);

    if any(missing_mask(:))
        fprintf('\nSecond pass: rerunning missing points with aggressive settings...\n');

        [ii_missing, jj_missing] = find(missing_mask);

        for kk = 1:numel(ii_missing)
            i = ii_missing(kk);
            j = jj_missing(kk);

            gamma = gamma_list(i);
            D     = D_list(j);

            fprintf('  Rerun missing point: D = %.4f, gamma = %.4f\n', D, gamma);

            opts_rerun = opts;

            if (D > opts.early_D_min) && (gamma < opts.early_gamma_max)
                % For the extreme corner, keep the short active-overlap window.
                opts_rerun.Xmax  = 500;
                opts_rerun.Nx    = 1000;
                opts_rerun.tEnd  = 30;
                opts_rerun.NtOut = 400;
                opts_rerun.gap   = 50;

                opts_rerun.overlap_tol = 1e-14;
                opts_rerun.MaxStep     = 0.15;
                opts_rerun.t_fit_start = 8;

                opts_rerun.local_centroid_half_window = 20;
                opts_rerun.interface_for_speed = 'localcentroid';
            else
                % For other missing points, use longer-time aggressive settings.
                opts_rerun.Xmax  = 500;
                opts_rerun.Nx    = 1000;
                opts_rerun.tEnd  = 450;
                opts_rerun.NtOut = 650;
                opts_rerun.gap   = 40;

                opts_rerun.overlap_tol = 1e-13;
                opts_rerun.MaxStep     = 0.4;

                opts_rerun.local_centroid_half_window = 20;
                opts_rerun.interface_for_speed = 'localcentroid';
            end

            try
                out2 = run_one_case_ode15s(D, gamma, opts_rerun);

                if isfinite(out2.c)
                    c_mat(i,j)     = out2.c;
                    slope_mat(i,j) = out2.slope;
                    R2_mat(i,j)    = out2.R2;
                    nFit_mat(i,j)  = out2.nFit;

                    fprintf('      success: c = %+ .6e, R2 = %.4f, nFit = %d, method = %s\n', ...
                        out2.c, out2.R2, out2.nFit, out2.fit_method);
                else
                    fprintf('      still missing with localcentroid; trying uvmax...\n');

                    opts_rerun.interface_for_speed = 'uvmax';

                    out3 = run_one_case_ode15s(D, gamma, opts_rerun);

                    if isfinite(out3.c)
                        c_mat(i,j)     = out3.c;
                        slope_mat(i,j) = out3.slope;
                        R2_mat(i,j)    = out3.R2;
                        nFit_mat(i,j)  = out3.nFit;

                        fprintf('      uvmax success: c = %+ .6e, R2 = %.4f, nFit = %d, method = %s\n', ...
                            out3.c, out3.R2, out3.nFit, out3.fit_method);
                    else
                        fprintf('      still missing after uvmax fallback.\n');
                    end
                end

            catch ME
                fprintf('      failed again: %s\n', ME.message);
            end
        end
    end

    %% ============================================================
    %  BASIC SUMMARY
    % ============================================================

    nFinite = nnz(isfinite(c_mat));
    fprintf('\nFinite c values found: %d / %d\n', nFinite, numel(c_mat));

    if nFinite == 0
        error(['No finite wave speeds were obtained. ', ...
               'Try reducing gap, increasing tEnd, or lowering overlap_tol.']);
    end

    fprintf('\nRange of raw c_mat:\n');
    fprintf('  min c = %.6f\n', min(c_mat(:), [], 'omitnan'));
    fprintf('  max c = %.6f\n', max(c_mat(:), [], 'omitnan'));
    fprintf('  number of NaN values = %d / %d\n', nnz(~isfinite(c_mat)), numel(c_mat));
    fprintf('  number of c > 3 values = %d / %d\n', nnz(c_mat(:) > 3), numel(c_mat));

    %% ============================================================
    %  MESHGRID
    % ============================================================

    [DD, GG] = meshgrid(D_list, gamma_list);

    %% ============================================================
    %  DIAGNOSTIC: PRINT MISSING POINTS
    % ============================================================

    missing_mask = ~isfinite(c_mat);

    if any(missing_mask(:))
        fprintf('\nMissing c values:\n');
        [ii, jj] = find(missing_mask);

        for kk = 1:numel(ii)
            fprintf('    D = %.4f, gamma = %.4f\n', ...
                D_list(jj(kk)), gamma_list(ii(kk)));
        end
    else
        fprintf('\nNo missing c values.\n');
    end

    %% ============================================================
    %  CONTOUR PLOT
    % ============================================================

    c_plot = fill_nan_for_plot(c_mat);

    dc_fill = 0.1;
    levels_fill = -0.5:dc_fill:3;

    figure;
    contourf(DD, GG, c_plot, levels_fill, 'LineColor', 'none');
    colorbar;
    clim([-0.5 3]);
    hold on;

    cmin_plot = min(c_plot(:), [], 'omitnan');
    cmax_plot = max(c_plot(:), [], 'omitnan');

    if cmin_plot <= 0 && cmax_plot >= 0
        contour(DD, GG, c_plot, [0 0], 'k', 'LineWidth', 2.0);
    end

    xlabel('$D$', 'Interpreter', 'latex');
    ylabel('$\gamma$', 'Interpreter', 'latex');
    title('Contour plot of signed wave speed $c(D,\gamma)$', 'Interpreter', 'latex');
    box on;
    grid off;

    %% ============================================================
    %  SURFACE PLOT
    % ============================================================

    figure;
    surf(DD, GG, c_mat, 'EdgeColor', 'none');
    colorbar;
    clim([-0.5 3]);

    xlabel('$D$', 'Interpreter', 'latex');
    ylabel('$\gamma$', 'Interpreter', 'latex');
    zlabel('$c$', 'Interpreter', 'latex');

    view(135, 30);

    xlim([min(D_list), max(D_list)]);
    ylim([0, max(gamma_list)]);
    zlim([-0.5 3]);

    box on;

    fprintf('\nDone.\n');

end

%% ========================================================================
%  SOLVE ONE PARAMETER CASE
% ========================================================================
function out = run_one_case_ode15s(D, gamma, opts)

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

    interface_for_speed = pick(opts, 'interface_for_speed', 'localcentroid');
    local_win = pick(opts, 'local_centroid_half_window', 12);

    RelTol  = pick(opts, 'RelTol', 1e-5);
    AbsTol  = pick(opts, 'AbsTol', 1e-7);
    MaxStep = pick(opts, 'MaxStep', 1.0);

    use_median_speed_fallback = pick(opts, 'use_median_speed_fallback', true);
    R2_fallback_threshold     = pick(opts, 'R2_fallback_threshold', 0.95);

    verbose = pick(opts, 'verbose', false);

    x  = linspace(0, Xmax, Nx).';
    dx = x(2) - x(1);

    Lmat = neumann_laplacian_1d(Nx, dx);

    xmid    = 0.5 * Xmax;
    x_uedge = xmid - gap/2;
    x_vedge = xmid + gap/2;

    [u0, v0] = initial_condition_gap(x, A, wIC, x_uedge, x_vedge);
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
    L_equal      = nan(Nt,1);

    prevL = NaN;

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

        Lk = find_interface_uv_equal(x, uk, vk, prevL, active_tol_equal);

        if ~isnan(Lk)
            L_equal(k) = Lk;
            prevL = Lk;
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
        case 'equal'
            Lfit = L_equal;
        otherwise
            error('Unknown interface_for_speed = %s', interface_for_speed);
    end

    fit = fit_interface_speed_active_window(t_sol, Lfit, overlap_max, overlap_tol, ...
                                            fit_start_eff, fit_last_fraction);

    fit_method = 'active-window-linear';

    if use_median_speed_fallback
        if ~isfinite(fit.c) || (~isfinite(fit.R2)) || fit.R2 < R2_fallback_threshold
            fit2 = fit_speed_by_median_local_velocity(t_sol, Lfit, fit_start_eff);

            if isfinite(fit2.c)
                fit = fit2;
                fit_method = 'median-local-velocity';
            end
        end
    end

    if verbose
        fprintf('D = %.4f, gamma = %.4f, c = %.6e, R2 = %.6f, method = %s\n', ...
            D, gamma, fit.c, fit.R2, fit_method);
    end

    out = struct();

    out.D = D;
    out.gamma = gamma;

    out.c = fit.c;
    out.slope = fit.slope;
    out.R2 = fit.R2;
    out.nFit = fit.nFit;
    out.fit_method = fit_method;

    out.t = t_sol;
    out.x = x;

    out.L_uvmax = L_uvmax;
    out.L_uvlocal = L_uvlocal;
    out.L_uvcentroid = L_uvcentroid;
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
%  ACTIVE-WINDOW FIT
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

%% ========================================================================
%  FIT ONE WINDOW
% ========================================================================
function fit = fit_one_window(t, L, idx)

    idx = idx(:);

    if numel(idx) < 2
        fit = empty_fit();
        return;
    end

    p = polyfit(t(idx), L(idx), 1);

    slope = p(1);
    intercept = p(2);
    c = -slope;

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

%% ========================================================================
%  FALLBACK: MEDIAN LOCAL SPEED
% ========================================================================
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

    c_local = -dL(good) ./ dt(good);
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
    fit.slope = -c_med;
    fit.intercept = NaN;
    fit.c = c_med;
    fit.R2 = NaN;
    fit.nFit = numel(c_late);
end

%% ========================================================================
%  EMPTY FIT
% ========================================================================
function fit = empty_fit()

    fit = struct();
    fit.slope = NaN;
    fit.intercept = NaN;
    fit.c = NaN;
    fit.R2 = NaN;
    fit.nFit = 0;
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

    if isempty(xcand)
        Lx = NaN;
        return;
    end

    if ~isnan(prevL)
        [~, ibest] = min(abs(xcand - prevL));
    else
        [~, ibest] = max(wcand);
    end

    Lx = xcand(ibest);
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

    p = polyfit(xx, ww, 2);

    if abs(p(1)) < 1e-14
        xpeak = x(idxMax);
        return;
    end

    xv = -p(2)/(2*p(1));

    if xv < xx(1) || xv > xx(3) || ~isfinite(xv)
        xpeak = x(idxMax);
    else
        xpeak = xv;
    end
end

%% ========================================================================
%  FILL NaNs FOR CONTOUR PLOT ONLY
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