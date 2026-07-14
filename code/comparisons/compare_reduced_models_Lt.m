function compare_reduced_models_Lt()
% compare_reduced_models_Lt
%
% Compare L(t) in the small- and large-diffusion regimes. The two L(t)
% comparisons use their own parameter sets, while a third figure collects
% the absolute errors from both regimes.

    clearvars -except ans;
    clc;
    close all;

    outDir = fileparts(mfilename('fullpath'));

    [small, large] = default_cases();

    fprintf_case_header('Small-diffusion regime', small);
    twoSmall = solve_two_species(small);
    maskSmall = comparison_mask(twoSmall, small);
    iSmall = find(maskSmall, 1, 'first');
    stefan = solve_reduced_stefan(small);
    LStefanRaw = interp1(stefan.t, stefan.L, twoSmall.t, 'linear', 'extrap');
    LStefan = align_at_time(twoSmall.t, LStefanRaw, ...
        twoSmall.t(iSmall), twoSmall.L(iSmall));

    fprintf_case_header('Large-diffusion regime', large);
    twoLarge = solve_two_species(large);
    maskLarge = comparison_mask(twoLarge, large);
    iLarge = find(maskLarge, 1, 'first');
    tAlignLarge = twoLarge.t(iLarge);
    LAlignLarge = twoLarge.L(iLarge);

    [cFlux, fluxInfo] = solve_corrected_flux_speed(large);
    LFlux = LAlignLarge + cFlux * (twoLarge.t - tAlignLarge);

    leadingSlope = sqrt(large.D) / (large.gamma_comp * sqrt(3));
    cFluxLeading = solve_speed_from_endpoint_slope( ...
        large.a, leadingSlope, large.c_min, large.c_max);
    LFluxLeading = LAlignLarge + cFluxLeading * (twoLarge.t - tAlignLarge);

    fprintf('\nLarge-D phase-plane speeds:\n');
    fprintf('  corrected flux: c = %.8f, V1_y(0) = %.8f\n', ...
        cFlux, fluxInfo.V1_y0);
    fprintf('  leading-order flux: c = %.8f\n', cFluxLeading);

    make_comparison_plots(twoSmall, LStefan, maskSmall, ...
        twoLarge, LFlux, LFluxLeading, maskLarge, outDir);
end

function [small, large] = default_cases()
    base.delta = 1e-4;
    base.a = 0.20;
    base.d_u = 0.5;
    base.d_v = 0.5;
    base.gap = 60;
    base.wIC = 0.01;

    base.Xmax = 120;
    base.Nx = 6000;
    base.tEnd = 160;
    base.NtOut = 801;

    base.RelTol = 1e-5;
    base.AbsTol = 1e-7;
    base.MaxStep = 0.25;
    base.interaction_tol = 1e-6;

    % Choose the late-time window after the initial transient.
    base.fit_t_min = 105;
    base.fit_t_max = 135;

    % One-species Stefan PDE settings.
    base.M_left = 120;
    base.Ny_reduced = 500;

    % Phase-plane speed search settings.
    base.c_min = -1;
    base.c_max = sqrt(2) * (0.5 - base.a) - 1e-5;
    base.largeD_bvp_ymax = 25;
    base.largeD_bvp_nmesh = 300;

    small = base;
    small.name = 'small_D_stefan';
    small.D = 0.01;
    small.gamma_comp = 1.0;
    small.kappa_stefan = small.gamma_comp;
    small.use_corrected_smallD = true;

    large = base;
    large.name = 'large_D_flux';
    large.D = 25.0;
    large.gamma_comp = 10.0;
    large.kappa_stefan = large.gamma_comp;
    large.use_corrected_smallD = false;
end

function fprintf_case_header(label, cfg)
    fprintf('\n============================================================\n');
    fprintf('%s\n', label);
    fprintf('D = %.6g, gamma = %.6g, a = %.6g, delta = %.3g\n', ...
        cfg.D, cfg.gamma_comp, cfg.a, cfg.delta);
    fprintf('============================================================\n');
end

function mask = comparison_mask(two, cfg)
    mask = two.t >= cfg.fit_t_min & two.t <= cfg.fit_t_max & isfinite(two.L);
    if nnz(mask) < 3
        validIdx = find(isfinite(two.L));
        if numel(validIdx) < 10
            error('Not enough valid L(t) points for case "%s".', cfg.name);
        end

        % If the interface leaves the finite computational domain before the
        % requested late-time window, use the latest available interior data.
        nSkipEnd = min(5, floor(0.05*numel(validIdx)));
        iEnd = numel(validIdx) - nSkipEnd;
        iStart = max(1, iEnd - 150);
        mask = false(size(two.t));
        mask(validIdx(iStart:iEnd)) = true;
        warning(['Case "%s" has no valid interface in the requested window. ' ...
            'Using the latest available interval t = [%.3g, %.3g].'], ...
            cfg.name, two.t(validIdx(iStart)), two.t(validIdx(iEnd)));
    end
end

function LAligned = align_at_time(t, L, tAlign, LAlign)
    offset = LAlign - interp1(t, L, tAlign, 'linear', 'extrap');
    LAligned = L + offset;
end

function sol = solve_two_species(cfg)
    x = linspace(0, cfg.Xmax, cfg.Nx).';
    dx = x(2) - x(1);
    Lmat = neumann_laplacian_1d(cfg.Nx, dx);

    xmid = 0.5 * cfg.Xmax;
    b_u = xmid - cfg.gap/2;
    b_v = xmid + cfg.gap/2;
    x_plot = x - b_u;

    [u0, v0] = initial_condition_gap(x, cfg.d_u, cfg.d_v, cfg.wIC, b_u, b_v);
    y0 = [u0; v0];
    tspan = linspace(0, cfg.tEnd, cfg.NtOut);

    rhs = @(t,y) two_species_rhs(t, y, Lmat, cfg.D, cfg.delta, cfg.gamma_comp, cfg.a);
    jac = @(t,y) two_species_jacobian(y, Lmat, cfg.D, cfg.delta, cfg.gamma_comp, cfg.a);
    Jpat = jacobian_pattern(cfg.Nx, Lmat);

    opts = odeset( ...
        'RelTol', cfg.RelTol, ...
        'AbsTol', cfg.AbsTol, ...
        'MaxStep', cfg.MaxStep, ...
        'Jacobian', jac, ...
        'JPattern', Jpat, ...
        'NonNegative', 1:(2*cfg.Nx));

    fprintf('Solving two-species PDE...\n');
    [t_sol, y_sol] = ode15s(rhs, tspan, y0, opts);
    fprintf('Two-species PDE done.\n');

    u = y_sol(:, 1:cfg.Nx);
    v = y_sol(:, cfg.Nx+1:end);
    uv_prod = u .* v;

    L_track = nan(numel(t_sol), 1);
    max_uv = nan(numel(t_sol), 1);
    for k = 1:numel(t_sol)
        [max_uv(k), idx] = max(uv_prod(k,:));
        if max_uv(k) > cfg.interaction_tol
            L_track(k) = subgrid_peak_location(x_plot, uv_prod(k,:), idx);
        end
    end

    sol.t = t_sol(:);
    sol.L = L_track(:);
    sol.max_uv = max_uv(:);
end

function sol = solve_reduced_stefan(cfg)
    y = linspace(0, 1, cfg.Ny_reduced).';
    dy = y(2) - y(1);

    L0 = 0;
    u0 = 0.5 * ones(cfg.Ny_reduced, 1);
    u0(end) = 0;

    % State vector stores u(1:end-1), with u(end)=0 imposed, and L.
    Y0 = [u0(1:end-1); L0];
    tspan = linspace(0, cfg.tEnd, cfg.NtOut);

    opts = odeset( ...
        'RelTol', 1e-5, ...
        'AbsTol', 1e-7, ...
        'MaxStep', cfg.MaxStep, ...
        'NonNegative', 1:(cfg.Ny_reduced-1));

    fprintf('Solving reduced Stefan one-species PDE...\n');
    [t_sol, Y_sol] = ode15s(@(t,Y) reduced_stefan_rhs(t, Y, cfg, dy, y), tspan, Y0, opts);
    fprintf('Reduced Stefan PDE done.\n');

    sol.t = t_sol(:);
    sol.L = Y_sol(:, end);
end

function dYdt = reduced_stefan_rhs(~, Y, cfg, dy, y)
    N = cfg.Ny_reduced;

    u = zeros(N, 1);
    u(1:N-1) = Y(1:N-1);
    u(N) = 0;

    L = Y(end);
    alpha = cfg.M_left + L;
    if alpha <= 1
        error('Reduced Stefan domain length became non-positive. Increase M_left.');
    end

    uy_right = (3*u(N) - 4*u(N-1) + u(N-2)) / (2*dy);
    q = -cfg.kappa_stefan / alpha * uy_right;

    if cfg.use_corrected_smallD
        % Corrected small-D condition:
        %   Ldot + D/(2 Ldot) = -gamma u_x(L^-).
        % Here q = -gamma u_x(L^-). The larger root recovers Ldot ~ q
        % as D -> 0. If the discriminant is negative during the early
        % transient, the corrected asymptotic relation is outside its
        % range of validity, so we use the limiting real value q/2.
        disc = q.^2 - 2*cfg.D;
        if disc >= 0
            Ldot = 0.5 * (q + sqrt(disc));
        else
            Ldot = 0.5 * q;
        end
    else
        Ldot = q;
    end

    uyy = zeros(N-1, 1);

    % Neumann condition at y=0: u_y=0 using a ghost point.
    uyy(1) = 2 * (u(2) - u(1)) / dy^2;

    for j = 2:N-2
        uyy(j) = (u(j+1) - 2*u(j) + u(j-1)) / dy^2;
    end

    % Last unknown adjacent to Dirichlet point u(N)=0.
    j = N - 1;
    uyy(j) = (u(j+1) - 2*u(j) + u(j-1)) / dy^2;

    uy = zeros(N-1, 1);
    uy(1) = 0;
    for j = 2:N-1
        uy(j) = (u(j+1) - u(j-1)) / (2*dy);
    end

    F = u(1:N-1) .* (1 - u(1:N-1)) .* (u(1:N-1) - cfg.a);
    dudt = (y(1:N-1) .* Ldot ./ alpha) .* uy + (1 / alpha^2) .* uyy + F;

    dYdt = [dudt; Ldot];
end

function c = solve_speed_from_endpoint_slope(a, targetSlope, cMin, cMax)
    % Find c such that -W(0;c) = targetSlope for the phase-plane equation.
    if targetSlope <= 0
        error('targetSlope must be positive.');
    end

    f = @(cc) endpoint_slope_magnitude(cc, a) - targetSlope;

    cGrid = linspace(cMin, cMax, 120);
    fGrid = nan(size(cGrid));
    for i = 1:numel(cGrid)
        try
            fGrid(i) = f(cGrid(i));
        catch
            fGrid(i) = NaN;
        end
    end

    idx = find(isfinite(fGrid(1:end-1)) & isfinite(fGrid(2:end)) & ...
        fGrid(1:end-1).*fGrid(2:end) <= 0, 1, 'last');

    if isempty(idx)
        [~, imin] = min(abs(fGrid));
        error(['Could not bracket reduced flux speed. Closest c = %.6g gives ' ...
            '-W(0) = %.6g, target = %.6g. Try changing gamma, D, or c_min.'], ...
            cGrid(imin), fGrid(imin) + targetSlope, targetSlope);
    end

    opts = optimset('TolX', 1e-8, 'Display', 'off');
    c = fzero(f, [cGrid(idx), cGrid(idx+1)], opts);
end

function [c, info] = solve_corrected_flux_speed(cfg)
    leadingSlope = sqrt(cfg.D) / (cfg.gamma_comp * sqrt(3));
    f = @(cc) corrected_flux_residual(cc, cfg, leadingSlope);

    cGrid = linspace(cfg.c_min, cfg.c_max, 80);
    fGrid = nan(size(cGrid));
    for i = 1:numel(cGrid)
        try
            fGrid(i) = f(cGrid(i));
        catch ME
            warning('Corrected flux residual failed at c = %.6g: %s', cGrid(i), ME.message);
            fGrid(i) = NaN;
        end
    end

    idx = find(isfinite(fGrid(1:end-1)) & isfinite(fGrid(2:end)) & ...
        fGrid(1:end-1).*fGrid(2:end) <= 0, 1, 'last');

    if isempty(idx)
        [~, imin] = min(abs(fGrid));
        [targetClosest, V1yClosest] = corrected_flux_target_slope(cGrid(imin), cfg, leadingSlope);
        error(['Could not bracket corrected reduced flux speed. Closest c = %.6g gives ' ...
            '-W(0) = %.6g, corrected target = %.6g, V1_y(0) = %.6g. ' ...
            'Try changing gamma, D, c_min, or c_max.'], ...
            cGrid(imin), fGrid(imin) + targetClosest, targetClosest, V1yClosest);
    end

    opts = optimset('TolX', 1e-8, 'Display', 'off');
    c = fzero(f, [cGrid(idx), cGrid(idx+1)], opts);

    [targetSlope, V1_y0] = corrected_flux_target_slope(c, cfg, leadingSlope);
    info.leading_slope = leadingSlope;
    info.V1_y0 = V1_y0;
    info.target_slope = targetSlope;
end

function r = corrected_flux_residual(c, cfg, leadingSlope)
    [targetSlope, ~] = corrected_flux_target_slope(c, cfg, leadingSlope);
    if targetSlope <= 0
        error('Corrected target slope must be positive, but got %.6g.', targetSlope);
    end
    r = endpoint_slope_magnitude(c, cfg.a) - targetSlope;
end

function [targetSlope, V1_y0] = corrected_flux_target_slope(c, cfg, leadingSlope)
    V1_y0 = solve_largeD_V1_y0(c, cfg.largeD_bvp_ymax, cfg.largeD_bvp_nmesh);

    % Corrected flux condition:
    %   u_x(L^-) = -sqrt(D)/(gamma*sqrt(3)) - V1_y(0)/gamma.
    % The phase-plane solver uses targetSlope = -u_x(L^-).
    targetSlope = leadingSlope + V1_y0 / cfg.gamma_comp;
end

function V1_y0 = solve_largeD_V1_y0(c, yMax, nMesh)
    % Solve the coupled truncated BVP for V0 and V1:
    %   V0'' + V0(1 - V0) = 0,
    %   V1'' + (1 - 2 V0) V1 = -c V0',
    % with V0(0)=0, V0(yMax)=1, V1(0)=0, V1(yMax)=0.
    y = linspace(0, yMax, nMesh);

    guess = @(yy) largeD_V1_initial_guess(yy);
    solinit = bvpinit(y, guess);
    opts = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NMax', max(5000, 10*nMesh));
    sol = bvp4c(@(yy, z) largeD_V1_ode(yy, z, c), @largeD_V1_bc, solinit, opts);

    z0 = deval(sol, 0);
    V1_y0 = z0(4);
end

function z = largeD_V1_initial_guess(y)
    V0 = tanh(y / sqrt(6));
    V0_y = (1 / sqrt(6)) * sech(y / sqrt(6)).^2;
    V1 = zeros(size(y));
    V1_y = zeros(size(y));
    z = [V0; V0_y; V1; V1_y];
end

function dzdy = largeD_V1_ode(~, z, c)
    V0 = z(1, :);
    V0_y = z(2, :);
    V1 = z(3, :);
    V1_y = z(4, :);

    S0 = V0 .* (1 - V0);
    Sp0 = 1 - 2*V0;

    dzdy = [ ...
        V0_y; ...
        -S0; ...
        V1_y; ...
        -Sp0 .* V1 - c .* V0_y];
end

function res = largeD_V1_bc(za, zb)
    res = [ ...
        za(1); ...
        zb(1) - 1; ...
        za(3); ...
        zb(3)];
end

function m = endpoint_slope_magnitude(c, a)
    F = @(U) U .* (1 - U) .* (U - a);
    Fp1 = a - 1;
    lambda = (-c + sqrt(c^2 - 4*Fp1)) / 2;

    eta = 1e-7;
    U0 = 1 - eta;
    V0 = -lambda * eta;

    Uspan = [U0, 1e-7];
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 2e-3);
    [~, V] = ode15s(@(U,V) (-c.*V - F(U)) ./ V, Uspan, V0, opts);

    m = -V(end);
end

function make_comparison_plots(twoSmall, LStefan, maskSmall, ...
        twoLarge, LFlux, LFluxLeading, maskLarge, outDir)
    validSmall = maskSmall & isfinite(twoSmall.L) & isfinite(LStefan);
    validLarge = maskLarge & isfinite(twoLarge.L) & isfinite(LFlux) & ...
        isfinite(LFluxLeading);

    colTwo = [0.80 0.20 0.10];
    colStefan = [0.00 0.20 0.85];
    colFlux = [0.00 0.20 0.85];
    colFluxLeading = [0.00 0.00 0.00];

    fig1 = figure('Color', 'w', 'Position', [100 120 760 500]);
    ax1 = axes(fig1);
    hold(ax1, 'on');

    hTwoSmall = plot(ax1, twoSmall.t(validSmall), twoSmall.L(validSmall), '-', ...
        'Color', colTwo, 'LineWidth', 2.0);
    hStefan = plot(ax1, twoSmall.t(validSmall), LStefan(validSmall), '--', ...
        'Color', colStefan, 'LineWidth', 1.9);

    xlabel(ax1, '$t$', 'Interpreter', 'latex', 'FontSize', 17);
    ylabel(ax1, '$L(t)$', 'Interpreter', 'latex', 'FontSize', 17);
    legend(ax1, [hTwoSmall, hStefan], ...
        {'Two-species model', 'Reduced Stefan-type model'}, ...
        'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
    grid(ax1, 'on');
    box(ax1, 'on');
    set(ax1, 'FontSize', 13, 'LineWidth', 1.0, ...
        'TickLabelInterpreter', 'latex', 'Color', 'w', ...
        'XColor', 'k', 'YColor', 'k');

    fig2 = figure('Color', 'w', 'Position', [900 120 760 500]);
    ax2 = axes(fig2);
    hold(ax2, 'on');

    hTwoLarge = plot(ax2, twoLarge.t(validLarge), twoLarge.L(validLarge), '-', ...
        'Color', colTwo, 'LineWidth', 2.0);
    hFlux = plot(ax2, twoLarge.t(validLarge), LFlux(validLarge), '--', ...
        'Color', colFlux, 'LineWidth', 1.9);
    hFluxLeading = plot(ax2, twoLarge.t(validLarge), ...
        LFluxLeading(validLarge), '--', ...
        'Color', colFluxLeading, 'LineWidth', 1.9);

    xlabel(ax2, '$t$', 'Interpreter', 'latex', 'FontSize', 17);
    ylabel(ax2, '$L(t)$', 'Interpreter', 'latex', 'FontSize', 17);
    legend(ax2, [hTwoLarge, hFlux, hFluxLeading], ...
        {'Two-species model', 'Reduced flux-type model', ...
         'Reduced leading-order flux-type model'}, ...
        'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
    grid(ax2, 'on');
    box(ax2, 'on');
    set(ax2, 'FontSize', 13, 'LineWidth', 1.0, ...
        'TickLabelInterpreter', 'latex', 'Color', 'w', ...
        'XColor', 'k', 'YColor', 'k');

    errorStefan = abs(LStefan - twoSmall.L);
    errorFlux = abs(LFlux - twoLarge.L);
    errorFluxLeading = abs(LFluxLeading - twoLarge.L);

    fig3 = figure('Color', 'w', 'Position', [500 650 760 500]);
    ax3 = axes(fig3);
    hold(ax3, 'on');

    validSmallErr = omit_first_valid_point(validSmall);
    validLargeErr = omit_first_valid_point(validLarge);

    hErrStefan = plot(ax3, twoSmall.t(validSmallErr), errorStefan(validSmallErr), '--', ...
        'Color', colStefan, 'LineWidth', 1.9);
    hErrFlux = plot(ax3, twoLarge.t(validLargeErr), errorFlux(validLargeErr), '-', ...
        'Color', [0.85 0.10 0.10], 'LineWidth', 1.9);
    hErrFluxLeading = plot(ax3, twoLarge.t(validLargeErr), ...
        errorFluxLeading(validLargeErr), '--', ...
        'Color', colFluxLeading, 'LineWidth', 1.9);

    xlabel(ax3, '$t$', 'Interpreter', 'latex', 'FontSize', 17);
    ylabel(ax3, 'Absolute errors', 'Interpreter', 'latex', 'FontSize', 17);
    legend(ax3, [hErrStefan, hErrFlux, hErrFluxLeading], ...
        {'Stefan-type error ($D=0.01$, $\gamma=1$)', ...
         'Corrected flux-type error ($D=25$, $\gamma=10$)', ...
         'Leading-order flux-type error ($D=25$, $\gamma=10$)'}, ...
        'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
    grid(ax3, 'on');
    box(ax3, 'on');
    set(ax3, 'FontSize', 13, 'LineWidth', 1.0, ...
        'TickLabelInterpreter', 'latex', 'Color', 'w', ...
        'XColor', 'k', 'YColor', 'k');

    savefig(fig1, fullfile(outDir, 'Lt_compare_smallD_stefan.fig'));
    savefig(fig2, fullfile(outDir, 'Lt_compare_largeD_flux.fig'));
    savefig(fig3, fullfile(outDir, 'Lt_error_reduced_models.fig'));
end

function maskOut = omit_first_valid_point(maskIn)
    maskOut = maskIn;
    firstIdx = find(maskOut, 1, 'first');
    if ~isempty(firstIdx) && nnz(maskOut) > 1
        maskOut(firstIdx) = false;
    end
end

function dydt = two_species_rhs(~, y, Lmat, D, delta, gamma_comp, a)
    Nx = numel(y) / 2;
    u = y(1:Nx);
    v = y(Nx+1:end);

    R = u .* (1-u) .* (u-a);
    S = v .* (1-v);

    dudt = Lmat*u + R - (u.*v)/delta;
    dvdt = D*(Lmat*v) + S - gamma_comp*(u.*v)/delta;

    dydt = [dudt; dvdt];
end

function J = two_species_jacobian(y, Lmat, D, delta, gamma_comp, a)
    Nx = numel(y) / 2;
    u = y(1:Nx);
    v = y(Nx+1:end);

    Rp = -3*u.^2 + 2*(1+a)*u - a;
    Sp = 1 - 2*v;

    J11 = Lmat + spdiags(Rp - v/delta, 0, Nx, Nx);
    J12 = spdiags(-u/delta, 0, Nx, Nx);
    J21 = spdiags(-gamma_comp*v/delta, 0, Nx, Nx);
    J22 = D*Lmat + spdiags(Sp - gamma_comp*u/delta, 0, Nx, Nx);

    J = [J11, J12; J21, J22];
end

function Jpat = jacobian_pattern(Nx, Lmat)
    S = spones(Lmat);
    I = speye(Nx);
    Jpat = [S I; I S];
end

function L = neumann_laplacian_1d(Nx, dx)
    e = ones(Nx, 1);
    L = spdiags([e -2*e e], -1:1, Nx, Nx) / dx^2;

    % Homogeneous Neumann boundary conditions via ghost points.
    L(1,1) = -2/dx^2;
    L(1,2) =  2/dx^2;
    L(Nx,Nx) = -2/dx^2;
    L(Nx,Nx-1) =  2/dx^2;
end

function [u0, v0] = initial_condition_gap(x, d_u, d_v, w, b_u, b_v)
    H_u = 0.5 * (1 + tanh((x - b_u) / w));
    H_v = 0.5 * (1 + tanh((x - b_v) / w));

    u0 = d_u * (1 - H_u);
    v0 = d_v * H_v;
end

function x_peak = subgrid_peak_location(x, f, idx)
    % Quadratic sub-grid interpolation of the local maximum of f.
    % This removes grid-point staircasing in L(t)=argmax_x u(x,t)v(x,t).
    x_peak = x(idx);

    if idx <= 1 || idx >= numel(f)
        return;
    end

    f_local = f([idx-1, idx, idx+1]);
    if any(~isfinite(f_local)) || any(f_local <= 0)
        return;
    end

    g = log(f_local);
    denom = g(1) - 2*g(2) + g(3);
    if abs(denom) < eps
        return;
    end

    offset = 0.5 * (g(1) - g(3)) / denom;
    offset = max(min(offset, 1), -1);
    dx = x(2) - x(1);
    x_peak = x(idx) + offset * dx;
end
