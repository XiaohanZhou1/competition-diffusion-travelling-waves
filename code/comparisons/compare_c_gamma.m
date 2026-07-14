function compare_c_gamma()
% compare_c_gamma
%
% Compare c(gamma) from the full two-species PDE with the corrected
% Stefan-type model for small D and the corrected/leading-order flux-type
% models for large D. A third figure collects the absolute speed errors.

    clearvars -except ans;
    clc;
    close all;

    a = 0.20;
    delta = 1e-4;

    D_stefan = 0.01;
    gamma_stefan = linspace(0.5, 1.6, 17);
    stefan = compute_stefan_curve(a, D_stefan);

    optsStefan = default_full_pde_options();
    optsStefan.delta = delta;
    optsStefan.a1 = a;
    optsStefan.interface_for_speed = 'ulevel';
    optsStefan.u_level = 0.05;
    Tstefan = compute_full_pde_curve(D_stefan, gamma_stefan, optsStefan, 'Stefan');

    D_flux = 25;
    eta0 = flux_zero_speed_eta(a);
    gamma_flux = linspace(0.55, 2.40, 31) * eta0 * sqrt(D_flux);
    flux = compute_flux_curves(a, D_flux);

    optsFlux = default_full_pde_options();
    optsFlux.delta = delta;
    optsFlux.a1 = a;
    optsFlux.A = 1.0;
    optsFlux.Xmax = 300;
    optsFlux.Nx = 650;
    optsFlux.tEnd = 220;
    optsFlux.NtOut = 360;
    optsFlux.gap = 80;
    optsFlux.wIC = 1.0;
    optsFlux.overlap_tol = 1e-9;
    optsFlux.t_fit_start = 20;
    optsFlux.fit_t_min = NaN;
    optsFlux.fit_t_max = NaN;
    optsFlux.interface_for_speed = 'localcentroid';
    optsFlux.RelTol = 1e-5;
    optsFlux.AbsTol = 1e-7;
    optsFlux.MaxStep = 0.8;
    optsFlux.use_median_speed_fallback = true;
    Tflux = compute_full_pde_curve(D_flux, gamma_flux, optsFlux, 'Flux');

    plot_comparisons(Tstefan, stefan.large, Tflux, flux, D_stefan, D_flux);
end

function Tfull = compute_full_pde_curve(D, gammaList, opts, label)
    rows = nan(numel(gammaList), 6);
    fprintf('\n%s full-PDE curve: D = %.6g, points = %d\n', ...
        label, D, numel(gammaList));

    for i = 1:numel(gammaList)
        fprintf('  case %d/%d: gamma = %.8f\n', i, numel(gammaList), gammaList(i));
        out = run_one_full_pde_case(D, gammaList(i), opts);
        rows(i,:) = [D, gammaList(i), out.c, out.slope, out.R2, out.nFit];
        fprintf('    c = %+.8e, R2 = %.4f, method = %s\n', ...
            out.c, out.R2, out.fit_method);
    end

    Tfull = array2table(rows, 'VariableNames', ...
        {'D','gamma','c','slopeL','R2','nFit'});
    Tfull = sortrows(Tfull, 'gamma');
end

function plot_comparisons(Tstefan, stefan, Tflux, flux, Dstefan, Dflux)
    red = [0.85 0.10 0.05];
    blue = [0.00 0.25 1.00];
    black = [0.00 0.00 0.00];

    [gammaStefanSmooth, cStefanFullSmooth] = smooth_full_curve(Tstefan);
    stefanMask = stefan.gamma >= min(Tstefan.gamma) & ...
        stefan.gamma <= max(Tstefan.gamma);

    fig1 = figure('Color', 'w', 'Position', [80 120 720 520]);
    ax1 = formatted_axes(fig1);
    hFull1 = plot(ax1, gammaStefanSmooth, cStefanFullSmooth, '-', ...
        'Color', red, 'LineWidth', 1.8);
    hReduced1 = plot(ax1, stefan.gamma(stefanMask), stefan.c(stefanMask), '--', ...
        'Color', blue, 'LineWidth', 2.0);
    finish_axes(ax1, [hFull1, hReduced1], ...
        {'Two-species model', 'Reduced Stefan-type model'});
    xlim(ax1, [min(Tstefan.gamma), max(Tstefan.gamma)]);

    [gammaFluxSmooth, cFluxFullSmooth] = smooth_full_curve(Tflux);
    fluxMask = flux.corrected.gamma >= min(Tflux.gamma) & ...
        flux.corrected.gamma <= max(Tflux.gamma);
    leadingMask = flux.leading.gamma >= min(Tflux.gamma) & ...
        flux.leading.gamma <= max(Tflux.gamma);

    fig2 = figure('Color', 'w', 'Position', [840 120 720 520]);
    ax2 = formatted_axes(fig2);
    hFull2 = plot(ax2, gammaFluxSmooth, cFluxFullSmooth, '-', ...
        'Color', red, 'LineWidth', 1.8);
    hReduced2 = plot(ax2, flux.corrected.gamma(fluxMask), ...
        flux.corrected.c(fluxMask), '--', 'Color', blue, 'LineWidth', 2.0);
    hLeading2 = plot(ax2, flux.leading.gamma(leadingMask), ...
        flux.leading.c(leadingMask), '--', 'Color', black, 'LineWidth', 1.8);
    finish_axes(ax2, [hFull2, hReduced2, hLeading2], ...
        {'Two-species model', 'Reduced flux-type model', ...
         'Reduced leading-order flux-type model'});
    xlim(ax2, [min(Tflux.gamma), max(Tflux.gamma)]);

    cStefanAtFull = interp1(stefan.gamma, stefan.c, Tstefan.gamma, 'pchip', NaN);
    cFluxAtFull = interp1(flux.corrected.gamma, flux.corrected.c, ...
        Tflux.gamma, 'pchip', NaN);
    cLeadingAtFull = interp1(flux.leading.gamma, flux.leading.c, ...
        Tflux.gamma, 'pchip', NaN);

    errStefan = abs(cStefanAtFull - Tstefan.c);
    errFlux = abs(cFluxAtFull - Tflux.c);
    errLeading = abs(cLeadingAtFull - Tflux.c);

    fig3 = figure('Color', 'w', 'Position', [460 680 720 520]);
    ax3 = formatted_axes(fig3);
    hErrStefan = plot(ax3, NaN, NaN, '-', ...
        'Color', blue, 'LineWidth', 1.8);
    hErrFlux = plot(ax3, Tflux.gamma, errFlux, '-', ...
        'Color', red, 'LineWidth', 1.8);
    hErrLeading = plot(ax3, Tflux.gamma, errLeading, '--', ...
        'Color', black, 'LineWidth', 1.8);
    xlabel(ax3, '$\gamma$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(ax3, '$|c_{\mathrm{reduced}}-c_{\mathrm{two}}|$', ...
        'Interpreter', 'latex', 'FontSize', 20);
    lgd = legend(ax3, [hErrStefan, hErrFlux, hErrLeading], ...
        {sprintf('Stefan-type error ($D=%.2g$)', Dstefan), ...
         sprintf('Corrected flux-type error ($D=%.2g$)', Dflux), ...
         sprintf('Leading-order flux-type error ($D=%.2g$)', Dflux)}, ...
        'Location', 'best', 'Interpreter', 'latex');
    format_legend(lgd);
    grid(ax3, 'on');
    box(ax3, 'on');

    insetMask = Tstefan.gamma >= 0.9 & Tstefan.gamma <= 1.5 & ...
        isfinite(errStefan);
    axInset = axes(fig3, 'Position', [0.20 0.56 0.30 0.30]);
    hold(axInset, 'on');
    plot(axInset, Tstefan.gamma(insetMask), errStefan(insetMask), '-', ...
        'Color', blue, 'LineWidth', 1.6);
    xlabel(axInset, '$\gamma$', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel(axInset, '$|\Delta c|$', 'Interpreter', 'latex', 'FontSize', 10);
    xlim(axInset, [0.9 1.5]);
    grid(axInset, 'on');
    box(axInset, 'on');
    set(axInset, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'FontSize', 9, 'TickLabelInterpreter', 'latex');
end

function ax = formatted_axes(fig)
    ax = axes(fig);
    hold(ax, 'on');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.65 0.65 0.65], 'MinorGridColor', [0.82 0.82 0.82], ...
        'FontSize', 12, 'TickLabelInterpreter', 'latex');
end

function finish_axes(ax, handles, labels)
    xlabel(ax, '$\gamma$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(ax, '$c$', 'Interpreter', 'latex', 'FontSize', 20);
    lgd = legend(ax, handles, labels, 'Location', 'best', 'Interpreter', 'latex');
    format_legend(lgd);
    grid(ax, 'on');
    box(ax, 'on');
end

function format_legend(lgd)
    lgd.FontSize = 12;
    lgd.Color = 'w';
    lgd.TextColor = 'k';
    lgd.EdgeColor = [0.2 0.2 0.2];
end

function [gammaSmooth, cSmooth] = smooth_full_curve(T)
    if height(T) >= 4
        gammaSmooth = linspace(min(T.gamma), max(T.gamma), 500).';
        cSmooth = pchip(T.gamma, T.c, gammaSmooth);
    else
        gammaSmooth = T.gamma;
        cSmooth = T.c;
    end
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
    opts.u_level = 0.05;
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
    u_level = pick(opts, 'u_level', 0.05);
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
    L_ulevel     = nan(Nt,1);

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

        L_ulevel(k) = descending_level_location(x, uk, u_level);
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
        case 'ulevel'
            Lfit = L_ulevel;
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
function stefan = compute_stefan_curve(a, D)
    F = @(U) U .* (1 - U) .* (U - a);
    cStar = sqrt(2) * (0.5 - a);

    cVals = [logspace(-6, -2, 120), ...
        linspace(1.05e-2, 0.98*cStar, 180), ...
        cStar - logspace(log10(0.02*cStar), -8, 120)].';
    cVals = unique(cVals, 'stable');
    cVals = cVals(cVals > 0 & cVals < cStar);

    kappaVals = nan(size(cVals));
    gammaVals = nan(size(cVals));
    for i = 1:numel(cVals)
        Wfront = front_slope_from_phase_plane(F, cVals(i));
        kappaVals(i) = -cVals(i) / Wfront;
        gammaVals(i) = -(cVals(i) + D/(2*cVals(i))) / Wfront;
    end

    keep = isfinite(gammaVals) & gammaVals > 0 & isfinite(cVals);
    gammaVals = gammaVals(keep);
    kappaVals = kappaVals(keep);
    cVals = cVals(keep);

    [cVals, idx] = sort(cVals);
    gammaVals = gammaVals(idx);
    kappaVals = kappaVals(idx);

    [~, iFold] = min(gammaVals);
    small = make_branch(gammaVals(1:iFold), cVals(1:iFold), kappaVals(1:iFold));
    large = make_branch(gammaVals(iFold:end), cVals(iFold:end), kappaVals(iFold:end));

    stefan = struct();
    stefan.small = small;
    stefan.large = large;
    stefan.gamma = [small.gamma; NaN; large.gamma];
    stefan.c = [small.c; NaN; large.c];
    stefan.kappa = [small.kappa; NaN; large.kappa];
end

function branch = make_branch(gammaVals, cVals, kappaVals)
    [gammaVals, idx] = sort(gammaVals);
    branch = struct();
    branch.gamma = gammaVals;
    branch.c = cVals(idx);
    branch.kappa = kappaVals(idx);
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

function flux = compute_flux_curves(a, D)
    F = @(U) U .* (1 - U) .* (U - a);
    cStar = sqrt(2) * (0.5 - a);
    yMax = 25;
    nMesh = 300;

    cVals = [linspace(1e-4, 0.98*cStar, 180), ...
        cStar - logspace(log10(0.02*cStar), -8, 120)].';
    cVals = unique(cVals, 'stable');
    cVals = cVals(cVals > 0 & cVals < cStar);

    gammaCorrected = nan(size(cVals));
    gammaLeading = nan(size(cVals));
    V1y0 = nan(size(cVals));

    for i = 1:numel(cVals)
        Wfront = front_slope_from_phase_plane(F, cVals(i));
        V1y0(i) = solve_largeD_V1_y0(cVals(i), yMax, nMesh);
        gammaCorrected(i) = -(sqrt(D)/sqrt(3) + V1y0(i)) / Wfront;
        gammaLeading(i) = -sqrt(D) / (sqrt(3) * Wfront);
    end

    flux.corrected = make_flux_branch(gammaCorrected, cVals);
    flux.leading = make_flux_branch(gammaLeading, cVals);
end

function branch = make_flux_branch(gammaVals, cVals)
    keep = isfinite(gammaVals) & gammaVals > 0 & isfinite(cVals);
    gammaVals = gammaVals(keep);
    cVals = cVals(keep);
    [gammaVals, idx] = sort(gammaVals);
    cVals = cVals(idx);
    [gammaVals, uniqueIdx] = unique(gammaVals, 'stable');
    branch.gamma = gammaVals;
    branch.c = cVals(uniqueIdx);
end

function eta0 = flux_zero_speed_eta(a)
    I0 = (1 - 2*a) / 12;
    eta0 = 1 / sqrt(6 * I0);
end

function V1_y0 = solve_largeD_V1_y0(c, yMax, nMesh)
    y = linspace(0, yMax, nMesh);
    solinit = bvpinit(y, @largeD_V1_initial_guess);
    opts = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
        'NMax', max(5000, 10*nMesh));
    sol = bvp4c(@(yy, z) largeD_V1_ode(yy, z, c), ...
        @largeD_V1_bc, solinit, opts);
    z0 = deval(sol, 0);
    V1_y0 = z0(4);
end

function z = largeD_V1_initial_guess(y)
    V0 = tanh(y / sqrt(6));
    V0_y = (1 / sqrt(6)) * sech(y / sqrt(6)).^2;
    z = [V0; V0_y; zeros(size(y)); zeros(size(y))];
end

function dzdy = largeD_V1_ode(~, z, c)
    V0 = z(1, :);
    V0_y = z(2, :);
    V1 = z(3, :);
    V1_y = z(4, :);
    dzdy = [V0_y; -V0 .* (1 - V0); V1_y; ...
        -(1 - 2*V0) .* V1 - c .* V0_y];
end

function res = largeD_V1_bc(za, zb)
    res = [za(1); zb(1) - 1; za(3); zb(3)];
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

function xlevel = descending_level_location(x, u, level)
    xlevel = NaN;
    u = u(:);
    x = x(:);

    idx = find(u(1:end-1) >= level & u(2:end) < level, 1, 'last');
    if isempty(idx)
        return;
    end

    u1 = u(idx);
    u2 = u(idx+1);
    if ~isfinite(u1) || ~isfinite(u2) || abs(u2-u1) < eps
        xlevel = x(idx);
        return;
    end

    theta = (level - u1) / (u2 - u1);
    theta = max(0, min(1, theta));
    xlevel = x(idx) + theta * (x(idx+1) - x(idx));
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
