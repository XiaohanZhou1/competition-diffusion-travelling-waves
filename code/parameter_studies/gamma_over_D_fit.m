function gamma_over_D_fit()
% smallD_phase_gamma_over_D_test
%
% Purpose:
%   Study the regime 0 < D <= 1 with a dense D-grid.
%   Output:
%     (1) a (D,gamma)-plane contourf plot of signed wave speed c(D,gamma)
%         with the zero-speed curve c=0 highlighted.
%     (2) a log-log power-law fit of the critical curve gamma_crit(D),
%         to test whether gamma_crit(D) ~ K*D^alpha, especially whether
%         alpha ~ 1 (equivalent to gamma/D ~ const).
%
% Convention:
%   L(t) ~ alpha0 + slope * t
%   c = -slope
%
% Hence:
%   c > 0  : forward travelling wave
%   c < 0  : backward travelling wave
%   c = 0  : critical curve

    clearvars -except ans;
    clc;
    close all;

    %% ============================================================
    %  USER PARAMETERS
    % ============================================================
    % Overall threshold scan with the original gap initial condition.
    % The main plot uses 0.1 < D < 15, with an inset for the small-D end.
    D_list = linspace(0.1, 15, 46);
    D_list = unique(D_list);

    gamma_list = linspace(0.02, 8, 25);

    % PDE / numerics
    opts = struct();
    opts.delta       = 1e-4;
    opts.a1          = 0.3;
    opts.A           = 1.0;

    opts.Xmax        = 250;
    opts.Nx          = 900;
    opts.tEnd        = 180;
    opts.Nt          = 360;

    % For very small D and short integrations, the original large gap can
    % prevent the two populations from interacting before tEnd.
    opts.gap         = 60;
    opts.wIC         = 0.01;
    opts.ic_type     = 'gap'; % 'gap' or 'stefan_small_branch'

    opts.eta_u       = 0.01;
    opts.eta_v       = 0.01;
    opts.active_tol  = 1e-4;

    opts.t_fit_start = 10;
    opts.verbose     = false;
    opts.interface_method = 'ulevel';
    opts.u_level = 0.05;
    opts.min_fit_points   = 12;
    opts.min_fit_span     = 8;

    opts.run_mode = 'overall_threshold_figures'; % 'overall_threshold_figures', 'small_speed_branch', or 'threshold'
    opts.min_forward_c = 1e-7;

    % Set false to use the original rectangular phase-diagram sweep.
    opts.use_bisection_only = true;
    opts.bracket_mode       = 'gamma_over_D';
    opts.gamma_over_D_bracket = [2, 150];
    opts.n_bracket_scan    = 5;
    opts.gamma_bracket      = [0.02, 8];
    opts.bisect_tol_c       = 1e-4;
    opts.bisect_tol_gamma   = 2e-4;
    opts.max_bisect_iter    = 20;
    opts.stop_on_c_tol      = true;
    opts.target_c           = 0.01;
    opts.output_dir         = fileparts(mfilename('fullpath'));

    %% ============================================================
    %  STORAGE
    % ============================================================
    ND     = numel(D_list);
    Ngamma = numel(gamma_list);

    if strcmpi(opts.run_mode, 'overall_threshold_figures')
        run_overall_threshold_figures(D_list, gamma_list, opts);
        return;
    end

    if strcmpi(opts.run_mode, 'small_speed_branch')
        run_small_speed_branch(D_list, gamma_list, opts);
        return;
    end

    if opts.use_bisection_only
        run_bisection_critical_curve(D_list, opts);
        return;
    end

    % rows = gamma, cols = D
    c_mat       = nan(Ngamma, ND);
    slope_mat   = nan(Ngamma, ND);
    R2_mat      = nan(Ngamma, ND);
    tc_mat      = nan(Ngamma, ND);
    fitStartMat = nan(Ngamma, ND);

    %% ============================================================
    %  PARAMETER SWEEP
    % ============================================================
    fprintf('Starting parameter sweep over small-D regime...\n');
    fprintf('Number of D values     = %d\n', ND);
    fprintf('Number of gamma values = %d\n', Ngamma);
    fprintf('Total PDE solves       = %d\n\n', ND * Ngamma);

    for j = 1:ND
        D = D_list(j);

        for i = 1:Ngamma
            gamma = gamma_list(i);

            fprintf('Case (%d/%d): D = %.4f, gamma = %.4f\n', ...
                (j-1)*Ngamma + i, ND*Ngamma, D, gamma);

            try
                out = run_one_case(D, gamma, opts);

                c_mat(i,j)       = out.c;
                slope_mat(i,j)   = out.slopeL;
                R2_mat(i,j)      = out.R2;
                tc_mat(i,j)      = out.tc;
                fitStartMat(i,j) = out.fit_start_eff;

            catch ME
                warning('Failed at D = %.6f, gamma = %.6f: %s', ...
                    D, gamma, ME.message);
            end
        end
    end

    fprintf('Finite speed values: %d / %d\n', nnz(isfinite(c_mat)), numel(c_mat));

    %% ============================================================
    %  GRID FOR PLOTTING
    % ============================================================
    [DD, GG] = meshgrid(D_list, gamma_list);

    % For contour plotting only, use a filled copy to reduce NaN breaks
    c_contour = c_mat;
    c_contour = fillmissing(c_contour, 'linear', 1, 'EndValues', 'none');
    c_contour = fillmissing(c_contour, 'linear', 2, 'EndValues', 'none');

    %% ============================================================
    %  FIGURE 1: (D,gamma)-PLANE
    %  Colour regions rich; 0 is explicitly a contour boundary
    % ============================================================
    figure;

    finite_c = c_contour(isfinite(c_contour));
    hold on;

    if numel(finite_c) >= 2 && min(finite_c) < max(finite_c)
        % Use manual fill levels so that c = 0 is exactly one colour boundary
        dc_fill = 0.05;
        cmin = min(finite_c);
        cmax = max(finite_c);

        low_level  = dc_fill * floor(cmin / dc_fill);
        high_level = dc_fill * ceil(cmax / dc_fill);
        levels_fill = low_level:dc_fill:high_level;

        if numel(levels_fill) < 2
            levels_fill = linspace(cmin, cmax, 12);
        end

        contourf(DD, GG, c_contour, levels_fill, 'LineColor', 'none');
        colorbar;

        % Highlight c = 0 only if the computed range crosses zero
        if cmin <= 0 && cmax >= 0
            contour(DD, GG, c_contour, [0 0], 'k', 'LineWidth', 2.2);
        end
    else
        warning('Not enough finite speed values to draw the contour plot.');
        plot(DD(isfinite(c_mat)), GG(isfinite(c_mat)), 'ko', 'MarkerFaceColor', 'k');
    end

    xlabel('D');
    ylabel('\gamma');
    grid off;
    box on;

    %% ============================================================
    %  ESTIMATE CRITICAL CURVE gamma_crit(D) FROM c = 0
    % ============================================================
    gamma_crit = nan(1, ND);

    for j = 1:ND
        c_col = c_mat(:,j);
        gamma_crit(j) = estimate_gamma_critical(gamma_list, c_col);
    end

    %% ============================================================
    %  FIGURE 2: POWER-LAW FIT OF gamma_crit(D)
    %  Test whether gamma_crit(D) ~ K * D^alpha, especially alpha ~ 1
    % ============================================================
    mask_valid = (D_list > 0) & ~isnan(gamma_crit) & (gamma_crit > 0);

    D_small = D_list(mask_valid);
    g_small = gamma_crit(mask_valid);

    if numel(D_small) >= 3
        % log-log least-squares fit
        p_small = polyfit(log(D_small(:)), log(g_small(:)), 1);
        alpha_small = p_small(1);
        K_small     = exp(p_small(2));

        gfit_small = K_small * D_small.^alpha_small;
        relerr_small = norm(g_small(:) - gfit_small(:)) / norm(g_small(:));

        fprintf('\n===== SMALL-D POWER-LAW FIT =====\n');
        fprintf('Fit over 0 < D <= 1\n');
        fprintf('gamma_crit(D) ~ K_small * D^(alpha_small)\n');
        fprintf('K_small     = %.8f\n', K_small);
        fprintf('alpha_small = %.8f\n', alpha_small);
        fprintf('relative L2 error = %.6e\n', relerr_small);

        % Direct test of gamma_crit(D)/D
        ratio_linear = g_small ./ D_small;
        fprintf('\n===== TEST OF gamma_crit(D)/D =====\n');
        fprintf('mean(gamma_crit/D) = %.8f\n', mean(ratio_linear, 'omitnan'));
        fprintf('std(gamma_crit/D)  = %.8f\n', std(ratio_linear, 'omitnan'));
        fprintf('min(gamma_crit/D)  = %.8f\n', min(ratio_linear));
        fprintf('max(gamma_crit/D)  = %.8f\n', max(ratio_linear));

        % smooth fit curve
        D_plot = logspace(log10(min(D_small)), log10(max(D_small)), 200);
        g_plot = K_small * D_plot.^alpha_small;

        figure;
        loglog(D_small, g_small, 'ko', 'LineWidth', 1.4, 'MarkerSize', 6, ...
            'DisplayName', '\gamma_{crit}(D)'); hold on;
        loglog(D_plot, g_plot, '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%.4f D^{%.4f}', K_small, alpha_small));

        xlabel('D');
        ylabel('\gamma_{crit}(D)');
        grid on;
        box on;

        lgd = legend('Location', 'best');
        lgd.ItemTokenSize = [18,10];
        lgd.FontSize = 11;

    else
        fprintf('\nNot enough valid critical points for the small-D power-law fit.\n');
    end

    %% ============================================================
    %  OPTIONAL DIAGNOSTIC: ratio gamma_crit(D)/D versus D
    %  (useful if you want a direct visual test of gamma/D ~ const)
    % ============================================================
    if numel(D_small) >= 3
        figure;
        plot(D_small, g_small ./ D_small, 'o-', 'LineWidth', 1.6, 'MarkerSize', 6);
        xlabel('D');
        ylabel('\gamma_{crit}(D)/D');
        grid on;
        box on;
    end

    %% ============================================================
    %  SUMMARY
    % ============================================================
    fprintf('\n=========== SUMMARY ===========\n');
    fprintf('Valid critical points found: %d / %d\n', nnz(mask_valid), ND);

    % Save tables if desired
    D_col        = DD(:);
    gamma_col    = GG(:);
    c_col        = c_mat(:);
    slope_col    = slope_mat(:);
    R2_col       = R2_mat(:);
    tc_col       = tc_mat(:);
    fitStart_col = fitStartMat(:);

    result_table = table(D_col, gamma_col, c_col, slope_col, R2_col, tc_col, fitStart_col, ...
        'VariableNames', {'D','gamma','c','slopeL','R2','tc','fitStart'});

    critical_table = table(D_list(:), gamma_crit(:), ...
        'VariableNames', {'D', 'gamma_crit'});

    disp(result_table(1:min(10,height(result_table)), :));
    disp(critical_table);

    % Uncomment if you want to save:
    % writetable(result_table,  'smallD_signed_c_results.csv');
    % writetable(critical_table,'smallD_critical_curve.csv');

    fprintf('Done.\n');
end

%% ========================================================================
%  OVERALL THRESHOLD FIGURES WITH ORIGINAL GAP INITIAL CONDITION
% ========================================================================
function run_overall_threshold_figures(D_list, ~, opts)
    opts.ic_type = 'gap';
    opts.interface_method = 'ulevel';

    ND = numel(D_list);
    gamma_crit = nan(size(D_list));
    c_crit = nan(size(D_list));

    fprintf('Starting fixed small-positive-speed scan with gap initial data...\n');
    fprintf('D range: %.4g to %.4g\n', min(D_list), max(D_list));
    fprintf('Number of D values = %d\n\n', ND);

    for j = 1:ND
        D = D_list(j);
        I0 = (1-2*opts.a1)/12;
        etaGuess = 1/sqrt(6*I0);
        gammaGuess = etaGuess*sqrt(D);
        bracket = [max(1e-4, 0.35*gammaGuess), 1.8*gammaGuess];

        fprintf('Threshold D = %.6g (%d/%d), bracket = [%.4g, %.4g]\n', ...
            D, j, ND, bracket(1), bracket(2));
        [gamma_crit(j), c_crit(j)] = find_gamma_by_bisection(D, bracket, opts);
        fprintf('  gamma_c0 ~= %.8g, c ~= %.4e\n\n', gamma_crit(j), c_crit(j));
    end

    targetC = pick(opts, 'target_c', 0);
    targetTol = 2*opts.bisect_tol_c;
    valid = isfinite(D_list) & isfinite(gamma_crit) & gamma_crit > 0 & ...
        isfinite(c_crit) & abs(c_crit-targetC) <= targetTol;
    Dv = D_list(valid);
    gv = gamma_crit(valid);

    largeFitMask = Dv > 0.1 & Dv < 15;
    if nnz(largeFitMask) >= 2
        xLarge = sqrt(Dv(largeFitMask));
        gLarge = gv(largeFitMask);
        pLargeLinear = polyfit(xLarge, gLarge, 1);
        slopeSqrtD = pLargeLinear(1);
        interceptSqrtD = pLargeLinear(2);
        relLarge = norm(gLarge - polyval(pLargeLinear, xLarge)) / norm(gLarge);

        fprintf('\n===== LARGE-D LINEAR FIT IN sqrt(D) =====\n');
        fprintf('gamma_crit(D) ~= %.8g sqrt(D) %+.8g\n', ...
            slopeSqrtD, interceptSqrtD);
        fprintf('fit range: 0.1 < D < 15\n');
        fprintf('relative L2 error = %.6e\n', relLarge);
    else
        slopeSqrtD = nan;
        interceptSqrtD = nan;
    end

    threshold_table = table(D_list(:), gamma_crit(:), c_crit(:), ...
        D_list(:)./sqrt(gamma_crit(:)), gamma_crit(:)./sqrt(D_list(:)), ...
        'VariableNames', {'D','gamma_crit','c_at_gamma_crit','D_over_sqrt_gamma','gamma_over_sqrtD'});

    %% Very-small-D threshold data for the inset linear fit
    D_small_fit = linspace(0.011, 0.019, 9);
    gamma_small_fit = nan(size(D_small_fit));
    c_small_fit = nan(size(D_small_fit));

    optsSmall = opts;
    optsSmall.ic_type = 'stefan_small_branch';
    optsSmall.interface_method = 'ulevel';
    optsSmall.u_level = 0.05;
    optsSmall.Xmax = 120;
    optsSmall.Nx = 6000;
    optsSmall.tEnd = 90;
    optsSmall.Nt = 451;
    optsSmall.t_fit_start = 5;
    optsSmall.min_fit_span = 10;

    fprintf('\nComputing Figure 2 small-D threshold fit...\n');
    for j = 1:numel(D_small_fit)
        D = D_small_fit(j);
        bracket = [0.02, 1.2];
        fprintf('Small-D threshold D = %.4g (%d/%d)\n', D, j, numel(D_small_fit));
        [gamma_small_fit(j), c_small_fit(j)] = ...
            find_gamma_by_bisection(D, bracket, optsSmall);
        fprintf('  gamma_c0 ~= %.8g, c ~= %.4e\n', ...
            gamma_small_fit(j), c_small_fit(j));
    end

    validSmall = isfinite(D_small_fit(:)) & isfinite(gamma_small_fit(:)) & ...
        (gamma_small_fit(:) > 0) & isfinite(c_small_fit(:)) & ...
        abs(c_small_fit(:)-targetC) <= targetTol;
    D_small_valid = D_small_fit(:);
    gamma_small_valid = gamma_small_fit(:);
    D_small_valid = D_small_valid(validSmall);
    gamma_small_valid = gamma_small_valid(validSmall);
    if numel(D_small_valid) >= 2
        pSmallLinear = polyfit(D_small_valid, gamma_small_valid, 1);
        slopeD = pSmallLinear(1);
        interceptSmall = pSmallLinear(2);
        relSmall = norm(gamma_small_valid - ...
            polyval(pSmallLinear, D_small_valid)) / norm(gamma_small_valid);

        fprintf('\n===== SMALL-D LINEAR FIT IN D =====\n');
        fprintf('gamma_c0(D) ~= %.8g D %+.8g\n', slopeD, interceptSmall);
        fprintf('fit range: 0.01 < D < 0.02\n');
        fprintf('relative L2 error = %.6e\n', relSmall);
        if abs(slopeD) > eps
            D10 = 10*abs(interceptSmall/slopeD);
            D05 = 20*abs(interceptSmall/slopeD);
            fprintf('empirical |intercept|/(|slope|D) < 10%% for D > %.8g\n', D10);
            fprintf('empirical |intercept|/(|slope|D) < 5%%  for D > %.8g\n', D05);
        end
    else
        slopeD = nan;
        interceptSmall = nan;
    end

    %% Full-PDE bridge points between the small- and main-D scans
    D_bridge_fit = [linspace(0.02, 0.05, 16), ...
        linspace(0.055, 0.095, 9)];
    gamma_bridge_fit = nan(size(D_bridge_fit));
    c_bridge_fit = nan(size(D_bridge_fit));

    optsBridge = optsSmall;
    fprintf('\nComputing full-PDE bridge points...\n');
    for j = 1:numel(D_bridge_fit)
        D = D_bridge_fit(j);
        bracket = [0.02, 1.5];
        fprintf('Bridge point D = %.4g (%d/%d)\n', ...
            D, j, numel(D_bridge_fit));
        [gamma_bridge_fit(j), c_bridge_fit(j)] = ...
            find_gamma_by_bisection(D, bracket, optsBridge);
        fprintf('  gamma_c0 ~= %.8g, c ~= %.4e\n', ...
            gamma_bridge_fit(j), c_bridge_fit(j));
    end

    validBridge = isfinite(D_bridge_fit(:)) & ...
        isfinite(gamma_bridge_fit(:)) & gamma_bridge_fit(:) > 0 & ...
        isfinite(c_bridge_fit(:)) & ...
        abs(c_bridge_fit(:)-targetC) <= targetTol;
    D_bridge_valid = D_bridge_fit(:);
    gamma_bridge_valid = gamma_bridge_fit(:);
    D_bridge_valid = D_bridge_valid(validBridge);
    gamma_bridge_valid = gamma_bridge_valid(validBridge);

    %% Both regressions shown on one common D-axis
    fig = figure('Color', 'w', 'Position', [120 120 780 540]);
    ax = axes(fig);
    hold(ax, 'on');

    Dtrue = [D_small_valid(:); D_bridge_valid(:); Dv(:)];
    gammaTrue = [gamma_small_valid(:); gamma_bridge_valid(:); gv(:)];
    [Dtrue, order] = sort(Dtrue);
    gammaTrue = gammaTrue(order);
    gammaForPlot = gammaTrue;
    bridgeSmoothMask = Dtrue >= 0.019 & Dtrue <= 0.05;
    if nnz(bridgeSmoothMask) >= 5
        gammaForPlot(bridgeSmoothMask) = smoothdata( ...
            gammaTrue(bridgeSmoothMask), 'sgolay', 5);
    end
    DtruePlot = linspace(min(Dtrue), max(Dtrue), 1000);
    gammaTruePlot = makima(Dtrue, gammaForPlot, DtruePlot);
    hNumerical = plot(ax, DtruePlot, gammaTruePlot, 'k-', ...
        'LineWidth', 1.7);

    DplotSmall = linspace(0.0101, 0.05, 300);
    hSmall = plot(ax, DplotSmall, ...
        slopeD*DplotSmall + interceptSmall, 'r--', 'LineWidth', 2.0);

    DplotLarge = linspace(0.02, 14.999, 500);
    hLarge = plot(ax, DplotLarge, ...
        slopeSqrtD*sqrt(DplotLarge) + interceptSqrtD, ...
        'b--', 'LineWidth', 2.2);

    xlabel(ax, '$D$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(ax, '$\gamma$', 'Interpreter', 'latex', 'FontSize', 16);
    legend(ax, [hNumerical, hSmall, hLarge], ...
        {'Numerical branch ($c=0.01$)', ...
         'Small-$D$ linear fit in $D$', ...
         'Large-$D$ linear fit in $\sqrt{D}$'}, ...
        'Location', 'best', 'Interpreter', 'latex');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
        'XScale', 'linear', 'YScale', 'linear', 'GridAlpha', 0.25);

    %% Small-D inset: the same curves restricted to D < 0.15, gamma < 0.9
    axInset = axes('Parent', fig, 'Position', [0.18 0.53 0.34 0.34]);
    hold(axInset, 'on');
    maskNumericalInset = DtruePlot <= 0.15;
    plot(axInset, DtruePlot(maskNumericalInset), ...
        gammaTruePlot(maskNumericalInset), 'k-', 'LineWidth', 1.4);
    maskSmallInset = DplotSmall <= 0.15;
    plot(axInset, DplotSmall(maskSmallInset), ...
        slopeD*DplotSmall(maskSmallInset) + interceptSmall, ...
        'r--', 'LineWidth', 1.5);
    maskLargeInset = DplotLarge <= 0.15;
    plot(axInset, DplotLarge(maskLargeInset), ...
        slopeSqrtD*sqrt(DplotLarge(maskLargeInset)) + interceptSqrtD, ...
        'b--', 'LineWidth', 1.5);
    xlim(axInset, [0 0.15]);
    ylim(axInset, [0 0.9]);
    xlabel(axInset, '$D$', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel(axInset, '$\gamma$', 'Interpreter', 'latex', 'FontSize', 10);
    grid(axInset, 'on');
    box(axInset, 'on');
    set(axInset, 'FontSize', 9, 'TickLabelInterpreter', 'latex', ...
        'GridAlpha', 0.22);

    %% Intercept-shifted full-range comparison
    validSmallPower = isfinite(D_small_valid) & isfinite(gamma_small_valid) & ...
        D_small_valid > 0 & gamma_small_valid > 0;
    validLargePower = isfinite(Dv) & isfinite(gv) & ...
        Dv > 0.1 & Dv < 15 & gv > 0;
    if nnz(validSmallPower) >= 2 && nnz(validLargePower) >= 2
        shiftedSmall = gamma_small_valid-interceptSmall;
        shiftedLarge = gv-interceptSmall;
        validShiftSmall = D_small_valid > 0 & shiftedSmall > 0 & ...
            isfinite(shiftedSmall);
        validShiftLarge = Dv > 0.1 & Dv < 15 & shiftedLarge > 0 & ...
            isfinite(shiftedLarge);

        if nnz(validShiftSmall) >= 2 && nnz(validShiftLarge) >= 2
            pShiftSmall = polyfit(log(D_small_valid(validShiftSmall)), ...
                log(shiftedSmall(validShiftSmall)),1);
            alphaShiftSmall = pShiftSmall(1);
            KShiftSmall = exp(pShiftSmall(2));
            pShiftLarge = polyfit(log(Dv(validShiftLarge)), ...
                log(shiftedLarge(validShiftLarge)),1);
            alphaShiftLarge = pShiftLarge(1);
            KShiftLarge = exp(pShiftLarge(2));

            shiftedForPlot = gammaForPlot-interceptSmall;
            validFullShift = Dtrue > 0 & shiftedForPlot > 0 & ...
                isfinite(shiftedForPlot);
            DshiftData = Dtrue(validFullShift);
            gammaShiftData = shiftedForPlot(validFullShift);
            DshiftPlot = logspace(log10(min(DshiftData)), ...
                log10(max(DshiftData)),1200);
            gammaShiftPlot = exp(makima(log(DshiftData), ...
                log(gammaShiftData),log(DshiftPlot)));

            DshiftSmallPlot = linspace(0.0101,0.07,300);
            DshiftLargePlot = linspace(0.05,14.999,500);

            %% Linear-axis representation of the shifted comparison
            figShift = figure('Color','w','Position',[170 130 800 550]);
            axShift = axes(figShift);
            hold(axShift,'on');
            hShiftNumerical = plot(axShift,DshiftPlot,gammaShiftPlot, ...
                'k-','LineWidth',1.8);
            hShiftSmall = plot(axShift,DshiftSmallPlot, ...
                KShiftSmall*DshiftSmallPlot.^alphaShiftSmall, ...
                'r--','LineWidth',2.0);
            hShiftLarge = plot(axShift,DshiftLargePlot, ...
                KShiftLarge*DshiftLargePlot.^alphaShiftLarge, ...
                'b--','LineWidth',2.1);
            xlabel(axShift,'$D$','Interpreter','latex','FontSize',16);
            ylabel(axShift,'$\gamma-b(0.01)$', ...
                'Interpreter','latex','FontSize',16);
            legend(axShift,[hShiftNumerical,hShiftSmall,hShiftLarge], ...
                {'Shifted numerical branch', ...
                sprintf('Small-$D$: $%.3gD^{%.4f}$', ...
                    KShiftSmall,alphaShiftSmall), ...
                sprintf('Large-$D$: $%.3gD^{%.4f}$', ...
                    KShiftLarge,alphaShiftLarge)}, ...
                'Location','best','Interpreter','latex');
            grid(axShift,'on'); box(axShift,'on');
            set(axShift,'FontSize',12,'TickLabelInterpreter','latex', ...
                'GridAlpha',0.25);

            %% Log-log representation of the same shifted comparison
            figShiftPower = figure('Color','w', ...
                'Position',[190 150 800 550]);
            axShiftPower = axes(figShiftPower);
            hold(axShiftPower,'on');
            hShiftNumericalPower = loglog(axShiftPower,DshiftPlot, ...
                gammaShiftPlot,'k-','LineWidth',1.8);
            hShiftSmallPower = loglog(axShiftPower,DshiftSmallPlot, ...
                KShiftSmall*DshiftSmallPlot.^alphaShiftSmall, ...
                'r--','LineWidth',2.0);
            hShiftLargePower = loglog(axShiftPower,DshiftLargePlot, ...
                KShiftLarge*DshiftLargePlot.^alphaShiftLarge, ...
                'b--','LineWidth',2.1);
            xlabel(axShiftPower,'$D$','Interpreter','latex','FontSize',16);
            ylabel(axShiftPower,'$\gamma-b(0.01)$', ...
                'Interpreter','latex','FontSize',16);
            legend(axShiftPower, ...
                [hShiftNumericalPower,hShiftSmallPower,hShiftLargePower], ...
                {'Shifted numerical branch', ...
                sprintf('Small-$D$: $%.3gD^{%.4f}$', ...
                    KShiftSmall,alphaShiftSmall), ...
                sprintf('Large-$D$: $%.3gD^{%.4f}$', ...
                    KShiftLarge,alphaShiftLarge)}, ...
                'Location','best','Interpreter','latex');
            grid(axShiftPower,'on'); box(axShiftPower,'on');
            set(axShiftPower,'FontSize',12, ...
                'TickLabelInterpreter','latex','XMinorGrid','on', ...
                'YMinorGrid','on','GridAlpha',0.25);

            fprintf('\n===== SHIFTED POWER-LAW CROSSOVER =====\n');
            fprintf('Shift b(c=0.01) = %.8g\n',interceptSmall);
            fprintf('Small-D shifted exponent = %.8g\n',alphaShiftSmall);
            fprintf('Large-D shifted exponent = %.8g\n',alphaShiftLarge);
        end
    end

    fprintf('\n=========== SUMMARY ===========\n');
    disp(threshold_table);
    if isfinite(slopeD)
        fprintf('Small-D fit: gamma_c0(D) ~= %.8g D %+.8g\n', ...
            slopeD, interceptSmall);
    end
    if isfinite(slopeSqrtD)
        fprintf('Large-D fit: gamma_crit(D) ~= %.8g sqrt(D) %+.8g\n', ...
            slopeSqrtD, interceptSqrtD);
    else
        fprintf('Not enough finite large-D threshold points for the D^(1/2) fit.\n');
    end
    disp(table(D_small_fit(:), gamma_small_fit(:), c_small_fit(:), ...
        'VariableNames', {'D','gamma_crit','c_at_gamma_crit'}));
    fprintf('\nFull-PDE bridge points:\n');
    disp(table(D_bridge_fit(:), gamma_bridge_fit(:), c_bridge_fit(:), ...
        'VariableNames', {'D','gamma_crit','c_at_gamma_crit'}));
    fprintf('Done.\n');
end

%% ========================================================================
%  SMALL-SPEED FORWARD BRANCH: c ~ K1*D/gamma
% ========================================================================
function run_small_speed_branch(D_list, gamma_list, opts)
    ND = numel(D_list);
    Ng = numel(gamma_list);

    c_mat       = nan(Ng, ND);
    c_branch_mat = nan(Ng, ND);
    slope_mat   = nan(Ng, ND);
    R2_mat      = nan(Ng, ND);
    tc_mat      = nan(Ng, ND);
    fitStartMat = nan(Ng, ND);

    fprintf('Starting small-speed forward branch scan...\n');
    fprintf('D values              = %d\n', ND);
    fprintf('gamma values          = %d\n', Ng);
    fprintf('Total PDE solves      = %d\n\n', ND * Ng);

    for j = 1:ND
        D = D_list(j);
        for i = 1:Ng
            gamma = gamma_list(i);
            fprintf('Case (%d/%d): D = %.4g, gamma = %.4g, D/gamma = %.4g\n', ...
                (j-1)*Ng + i, ND*Ng, D, gamma, D/gamma);

            try
                out = run_one_case(D, gamma, opts);
                c_mat(i,j)       = out.c;
                c_branch_mat(i,j) = out.slopeL;
                slope_mat(i,j)   = out.slopeL;
                R2_mat(i,j)      = out.R2;
                tc_mat(i,j)      = out.tc;
                fitStartMat(i,j) = out.fit_start_eff;
            catch ME
                warning('Failed at D = %.6g, gamma = %.6g: %s', D, gamma, ME.message);
            end
        end
    end

    [DD, GG] = meshgrid(D_list, gamma_list);
    x_col = DD(:) ./ GG(:);
    c_col = c_branch_mat(:);
    forward = isfinite(c_col) & c_col > pick(opts, 'min_forward_c', 0);

    if nnz(forward) >= 2
        p0 = sum(x_col(forward).*c_col(forward)) / sum(x_col(forward).^2);
        p1 = polyfit(x_col(forward), c_col(forward), 1);
        c_fit0 = p0 * x_col(forward);
        relerr0 = norm(c_col(forward) - c_fit0) / norm(c_col(forward));

        fprintf('\n===== SMALL-SPEED BRANCH FIT =====\n');
        fprintf('Fit target: c_branch = L_t ~= K1*(D/gamma), using forward finite points only.\n');
        fprintf('Number of fit points = %d / %d\n', nnz(forward), numel(c_col));
        fprintf('Zero-intercept K1    = %.8g\n', p0);
        fprintf('Zero-intercept relative L2 error = %.6e\n', relerr0);
        fprintf('Free fit: c_branch ~= %.8g*(D/gamma) + %.8g\n', p1(1), p1(2));
    else
        p0 = nan;
        p1 = [nan, nan];
        relerr0 = nan;
        fprintf('\nNot enough forward finite points to fit c_branch ~= K1*(D/gamma).\n');
    end

    rows = table(DD(:), GG(:), x_col, c_mat(:), c_col, slope_mat(:), R2_mat(:), ...
        tc_mat(:), fitStartMat(:), ...
        'VariableNames', {'D','gamma','D_over_gamma','signed_c_old','c_branch_Lt','slopeL','R2','tc','fitStart'});
    rows.forward_fit_point = forward;

    figure('Color','w');
    hold on;
    for i = 1:Ng
        D_over_gamma = D_list ./ gamma_list(i);
        y = c_branch_mat(i,:);
        valid = isfinite(y);
        plot(D_over_gamma(valid), y(valid), 'o-', 'LineWidth', 1.5, ...
            'MarkerSize', 5, 'DisplayName', sprintf('\\gamma = %.3g', gamma_list(i)));
    end
    if isfinite(p0)
        xfit = linspace(0, max(x_col(forward)), 200);
        plot(xfit, p0*xfit, 'k--', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('fit: c = %.3g D/\\gamma', p0));
    end
    xlabel('D/\gamma');
    ylabel('c');
    grid on; box on;
    legend('Location','best');

    figure('Color','w');
    hold on;
    for i = 1:Ng
        y = gamma_list(i) * c_branch_mat(i,:);
        valid = isfinite(y);
        plot(D_list(valid), y(valid), 'o-', 'LineWidth', 1.5, ...
            'MarkerSize', 5, 'DisplayName', sprintf('\\gamma = %.3g', gamma_list(i)));
    end
    if isfinite(p0)
        xfit = linspace(0, max(DD(forward)), 200);
        plot(xfit, p0*xfit, 'k--', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('fit: \\gamma c = %.3g D', p0));
    end
    xlabel('D');
    ylabel('\gamma c');
    grid on; box on;
    legend('Location','best');

    fprintf('\n=========== SUMMARY ===========\n');
    disp(rows);
    fprintf('Done.\n');

end

%% ========================================================================
%  BISECTION-BASED CRITICAL CURVE
% ========================================================================
function run_bisection_critical_curve(D_list, opts)
    gamma_crit = nan(size(D_list));
    c_crit = nan(size(D_list));

    fprintf('Starting bisection scan for very small D.\n');
    fprintf('Number of D values = %d\n', numel(D_list));
    fprintf('Each D uses at most %d bisection iterations.\n\n', opts.max_bisect_iter);

    gammaGuess = nan;
    for j = 1:numel(D_list)
        D = D_list(j);

        if strcmpi(pick(opts, 'bracket_mode', 'absolute'), 'gamma_over_D')
            ratioBracket = pick(opts, 'gamma_over_D_bracket', [20, 120]);
            bracket = D * ratioBracket;
        elseif isfinite(gammaGuess)
            bracket = [max(1e-4, 0.65*gammaGuess), 1.45*gammaGuess];
        else
            bracket = opts.gamma_bracket;
        end

        fprintf('D = %.6g (%d/%d), bracket = [%.5g, %.5g]\n', ...
            D, j, numel(D_list), bracket(1), bracket(2));

        [gamma_crit(j), c_crit(j)] = find_gamma_by_bisection(D, bracket, opts);

        if isfinite(gamma_crit(j)) && ~strcmpi(pick(opts, 'bracket_mode', 'absolute'), 'gamma_over_D')
            gammaGuess = gamma_crit(j);
        end

        fprintf('  gamma_crit ~= %.8g, c ~= %.4e\n\n', gamma_crit(j), c_crit(j));
    end

    valid = isfinite(D_list) & isfinite(gamma_crit) & gamma_crit > 0;
    D = D_list(valid);
    g = gamma_crit(valid);

    if numel(D) >= 3
        fitSqrt = polyfit(log(D(:)), log(g(:)), 1);
        alpha = fitSqrt(1);
        K = exp(fitSqrt(2));

        fprintf('Power-law fit over bisection data:\n');
        fprintf('  gamma_crit(D) ~= %.8g D^{%.8g}\n', K, alpha);
        fprintf('  mean gamma/D = %.8g\n', mean(g ./ D));
        fprintf('  mean gamma/sqrt(D) = %.8g\n', mean(g ./ sqrt(D)));
    end

    figure('Color','w');
    loglog(D, g, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); hold on;
    if numel(D) >= 3
        Dplot = logspace(log10(min(D)), log10(max(D)), 200);
        loglog(Dplot, K * Dplot.^alpha, '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('Fit: %.4g D^{%.4f}', K, alpha));
    end
    xlabel('$D$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\gamma_{\rm crit}(D)$', 'Interpreter', 'latex', 'FontSize', 16);
    grid on; box on;
    legend({'$\gamma_{\rm crit}(D)$ data', 'Power-law fit'}, ...
        'Interpreter', 'latex', 'Location', 'best');
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

    figure('Color','w');
    plot(D, g ./ D, 'o-', 'LineWidth', 1.6, 'MarkerSize', 6); hold on;
    xlabel('$D$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\gamma_{\rm crit}(D)/D$', 'Interpreter', 'latex', 'FontSize', 16);
    grid on; box on;
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

    figure('Color','w');
    plot(D, g ./ sqrt(D), 'o-', 'LineWidth', 1.6, 'MarkerSize', 6); hold on;
    xlabel('$D$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\gamma_{\rm crit}(D)/\sqrt{D}$', 'Interpreter', 'latex', 'FontSize', 16);
    grid on; box on;
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

    result_table = table(D(:), g(:), c_crit(valid).', g(:)./D(:), g(:)./sqrt(D(:)), ...
        'VariableNames', {'D','gamma_crit','c_at_gamma_crit','gamma_over_D','gamma_over_sqrtD'});
    disp(result_table);

end

function [gammaRoot, cRoot] = find_gamma_by_bisection(D, bracket, opts)
    targetC = pick(opts, 'target_c', 0);
    [bracket, scanTable] = find_finite_sign_bracket(D, bracket, opts);
    if ~isempty(scanTable)
        fprintf('  coarse scan finite values:\n');
        disp(scanTable);
    end

    if any(~isfinite(bracket))
        warning('Could not find a finite sign-change bracket for D=%.6g.', D);
        gammaRoot = nan;
        cRoot = nan;
        return;
    end

    gL = bracket(1);
    gR = bracket(2);
    cL = run_one_case(D, gL, opts).c;
    cR = run_one_case(D, gR, opts).c;
    rL = cL - targetC;
    rR = cR - targetC;
    fprintf('  endpoint: gamma = %.7g, c = %.4e\n', gL, cL);
    fprintf('  endpoint: gamma = %.7g, c = %.4e\n', gR, cR);
    [bestGamma, bestC] = best_finite_speed([gL, gR], [cL, cR], targetC);

    expandCount = 0;
    while (~isfinite(rL) || ~isfinite(rR) || rL*rR > 0) && expandCount < 8
        if ~isfinite(cL)
            gL = max(1e-8, 0.5*gL);
            cL = run_one_case(D, gL, opts).c;
            rL = cL - targetC;
            fprintf('  expand L: gamma = %.7g, c = %.4e\n', gL, cL);
        elseif ~isfinite(cR)
            gR = 1.5*gR;
            cR = run_one_case(D, gR, opts).c;
            rR = cR - targetC;
            fprintf('  expand R: gamma = %.7g, c = %.4e\n', gR, cR);
        elseif rL < 0 && rR < 0
            gL = max(1e-4, 0.5*gL);
            cL = run_one_case(D, gL, opts).c;
            rL = cL - targetC;
            fprintf('  expand L: gamma = %.7g, c = %.4e\n', gL, cL);
        else
            gR = 1.5*gR;
            cR = run_one_case(D, gR, opts).c;
            rR = cR - targetC;
            fprintf('  expand R: gamma = %.7g, c = %.4e\n', gR, cR);
        end
        expandCount = expandCount + 1;
    end

    if ~isfinite(rL) || ~isfinite(rR) || rL*rR > 0
        warning('Could not bracket c=%.6g for D=%.6g.', targetC, D);
        gammaRoot = nan;
        cRoot = nan;
        return;
    end

    for iter = 1:opts.max_bisect_iter
        gM = 0.5*(gL + gR);
        cM = run_one_case(D, gM, opts).c;
        rM = cM - targetC;
        fprintf('    iter %02d: gamma = %.7g, c = %.4e\n', iter, gM, cM);

        if ~isfinite(cM)
            gL = gM;
            continue;
        end

        if ~isfinite(bestC) || abs(rM) < abs(bestC-targetC)
            bestGamma = gM;
            bestC = cM;
        end

        stopOnC = pick(opts, 'stop_on_c_tol', true);
        if (stopOnC && abs(rM) < opts.bisect_tol_c) || abs(gR - gL) < opts.bisect_tol_gamma
            gammaRoot = gM;
            cRoot = cM;
            return;
        end

        if rL*rM <= 0
            gR = gM;
            cR = cM; %#ok<NASGU>
            rR = rM; %#ok<NASGU>
        else
            gL = gM;
            cL = cM;
            rL = rM;
        end
    end

    gammaRoot = 0.5*(gL + gR);
    cRoot = run_one_case(D, gammaRoot, opts).c;
    if isfinite(bestGamma) && ...
            (~isfinite(cRoot) || abs(bestC-targetC) < abs(cRoot-targetC))
        gammaRoot = bestGamma;
        cRoot = bestC;
    end
end

function [bestGamma, bestC] = best_finite_speed(g, c, targetC)
    finite = isfinite(c);
    if ~any(finite)
        bestGamma = nan;
        bestC = nan;
        return;
    end

    gf = g(finite);
    cf = c(finite);
    [~, idx] = min(abs(cf-targetC));
    bestGamma = gf(idx);
    bestC = cf(idx);
end

%% ========================================================================
%  SOLVE ONE PARAMETER CASE
% ========================================================================
function out = run_one_case(D, gamma, opts)

    % ===== Defaults =====
    delta       = pick(opts, 'delta',       1e-3);
    a1          = pick(opts, 'a1',          0.1);
    A           = pick(opts, 'A',           1.0);

    Xmax        = pick(opts, 'Xmax',        250);
    Nx          = pick(opts, 'Nx',          800);
    tEnd        = pick(opts, 'tEnd',        160);
    Nt          = pick(opts, 'Nt',          360);

    gap         = pick(opts, 'gap',         80);
    wIC         = pick(opts, 'wIC',         1.0);
    ic_type     = pick(opts, 'ic_type',     'gap');

    eta_u       = pick(opts, 'eta_u',       0.01);
    eta_v       = pick(opts, 'eta_v',       0.01);
    active_tol  = pick(opts, 'active_tol',  1e-6);

    t_fit_start = pick(opts, 't_fit_start', 10);
    verbose     = pick(opts, 'verbose',     false);
    interface_method = pick(opts, 'interface_method', 'reaction_peak');
    local_centroid_half_window = pick(opts, 'local_centroid_half_window', 12);
    u_level = pick(opts, 'u_level', 0.05);
    min_fit_points   = pick(opts, 'min_fit_points', 8);
    min_fit_span     = pick(opts, 'min_fit_span', 0);

    % ===== Domain/time =====
    x    = linspace(0, Xmax, Nx).';
    dx   = x(2) - x(1);
    t    = linspace(0, tEnd, Nt);

    % ===== Initial conditions =====
    xmid    = 0.5 * Xmax;
    x_uedge = xmid - gap/2;
    x_vedge = xmid + gap/2;

    switch lower(ic_type)
        case 'gap'
            [u0, v0] = initial_condition_gap(x, A, wIC, x_uedge, x_vedge);
        case 'stefan_small_branch'
            L0 = xmid;
            [u0, v0] = initial_condition_stefan_small_branch(x, A, wIC, L0, a1, D);
        otherwise
            error('Unknown ic_type: %s', ic_type);
    end
    y0 = [u0; v0];

    % ===== Solve PDE using method of lines + ode15s =====
    Lmat = neumann_laplacian_1d(Nx, dx);
    rhs = @(t,y) rhs_fun(t, y, Lmat, D, delta, gamma, a1);
    jac = @(t,y) jacobian_fun(y, Lmat, D, delta, gamma, a1);
    Jpat = jacobian_pattern(Nx, Lmat);

    ode_opts = odeset( ...
        'RelTol', 1e-5, ...
        'AbsTol', 1e-7, ...
        'MaxStep', 0.25, ...
        'Jacobian', jac, ...
        'JPattern', Jpat, ...
        'NonNegative', 1:(2*Nx));

    [tcol, y_sol] = ode15s(rhs, t, y0, ode_opts);
    u = y_sol(:, 1:Nx);
    v = y_sol(:, Nx+1:2*Nx);
    x = x.';

    % ===== Track Xu(t), Xv(t) =====
    nTime = numel(tcol);
    Xu = nan(nTime,1);
    Xv = nan(nTime,1);

    for k = 1:nTime
        idx_u = find(u(k,:) >= eta_u, 1, 'last');
        if ~isempty(idx_u)
            Xu(k) = x(idx_u);
        end

        idx_v = find(v(k,:) >= eta_v, 1, 'first');
        if ~isempty(idx_v)
            Xv(k) = x(idx_v);
        end
    end

    % ===== Detect contact time =====
    tc_idx = find(Xu >= Xv, 1, 'first');
    if strcmpi(interface_method, 'ulevel')
        tc_idx = 1;
    end

    if isempty(tc_idx)
        tc = NaN;
    else
        tc = t(tc_idx);
    end

    % ===== Define interface L(t) by u=v =====
    L = nan(nTime,1);

    if ~isempty(tc_idx)
        prevL = NaN;
        for k = tc_idx:nTime
            switch lower(interface_method)
                case 'reaction_peak'
                    Lk = find_interface_reaction_peak(x, u(k,:), v(k,:), prevL, active_tol);
                case 'localcentroid'
                    Lk = find_interface_local_centroid(x, u(k,:), v(k,:), ...
                        prevL, active_tol^2, local_centroid_half_window);
                case 'ulevel'
                    Lk = descending_level_location(x, u(k,:), u_level);
                case 'uv_equal'
                    Lk = find_interface_uv_equal(x, u(k,:), v(k,:), prevL, active_tol);
                otherwise
                    error('Unknown interface_method: %s', interface_method);
            end
            L(k) = Lk;
            if ~isnan(Lk)
                prevL = Lk;
            end
        end
    end

    % ===== Fit linear phase =====
    if isnan(tc)
        fit_start_eff = t_fit_start;
    else
        fit_start_eff = max(t_fit_start, tc);
    end

    valid_idx = find((tcol > fit_start_eff) & ~isnan(L));

    if isempty(valid_idx)
        valid_span = 0;
    else
        valid_span = max(tcol(valid_idx)) - min(tcol(valid_idx));
    end

    if numel(valid_idx) >= min_fit_points && valid_span >= min_fit_span
        % Fit only the last 40% of valid points
        nValid = numel(valid_idx);
        iStart = max(1, floor(0.6*nValid));
        fit_idx = valid_idx(iStart:end);

        fit_mask = false(size(tcol));
        fit_mask(fit_idx) = true;

        p = polyfit(tcol(fit_idx), L(fit_idx), 1);

        slopeL     = p(1);
        interceptL = p(2);

        % Signed wave speed: c = dL/dt.
        c = slopeL;

        L_fit = polyval(p, tcol(fit_idx));
        ydata = L(fit_idx);

        SSres = sum((ydata - L_fit).^2);
        SStot = sum((ydata - mean(ydata)).^2);

        if SStot > 0
            R2 = 1 - SSres/SStot;
        else
            R2 = NaN;
        end

    elseif numel(valid_idx) >= max(2, min_fit_points) && valid_span >= min_fit_span
        fit_idx = valid_idx;

        fit_mask = false(size(tcol));
        fit_mask(fit_idx) = true;

        p = polyfit(tcol(fit_idx), L(fit_idx), 1);

        slopeL     = p(1);
        interceptL = p(2);
        c          = slopeL;

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
        fit_mask   = false(size(tcol));
        slopeL     = NaN;
        interceptL = NaN;
        c          = NaN;
        R2         = NaN;
    end

    if verbose
        fprintf('D = %.4f, gamma = %.4f, c = %.6f, R2 = %.6f\n', ...
            D, gamma, c, R2);
    end

    out = struct();
    out.D             = D;
    out.gamma         = gamma;
    out.tc            = tc;
    out.x             = x;
    out.t             = tcol;
    out.Xu            = Xu;
    out.Xv            = Xv;
    out.L             = L;
    out.slopeL        = slopeL;
    out.interceptL    = interceptL;
    out.c             = c;
    out.R2            = R2;
    out.fit_mask      = fit_mask;
    out.fit_start_eff = fit_start_eff;
end

%% ========================================================================
%  ESTIMATE CRITICAL gamma FOR ONE FIXED D
% ========================================================================
function gamma_crit = estimate_gamma_critical(gamma_list, c_col)
% For one fixed D, estimate gamma_crit where c crosses zero.
% Linear interpolation between neighboring gamma values with opposite signs.

    gamma_crit = NaN;

    valid = ~isnan(c_col);
    if nnz(valid) < 2
        return;
    end

    g = gamma_list(valid);
    c = c_col(valid);

    % exact zero
    [cmin, idxmin] = min(abs(c));
    if cmin < 1e-10
        gamma_crit = g(idxmin);
        return;
    end

    % sign-change intervals
    idx = find(c(1:end-1).*c(2:end) < 0);

    if isempty(idx)
        return;
    end

    % choose the sign-change interval nearest to zero in local magnitude
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

    gamma_crit = g1 - c1 * (g2 - g1) / (c2 - c1);
end

%% ========================================================================
%  LOCAL UTILITIES
% ========================================================================
function val = pick(s, name, defaultVal)
    if isfield(s, name)
        val = s.(name);
    else
        val = defaultVal;
    end
end

function [K, alpha] = fit_power_law(D, gamma)
    valid = isfinite(D) & isfinite(gamma) & D > 0 & gamma > 0;
    if nnz(valid) < 2
        K = nan;
        alpha = nan;
        return;
    end

    p = polyfit(log(D(valid)), log(gamma(valid)), 1);
    alpha = p(1);
    K = exp(p(2));
end

function [c,f,s] = pdefun(x,t,U,DUdx,D,delta,gamma,R,S) %#ok<INUSD>
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
    H1 = 0.5 * (1 + tanh((x - x_uedge)/w));
    H2 = 0.5 * (1 + tanh((x - x_vedge)/w));

    u0 = A * (1 - H1);
    v0 = A * H2;

    U0 = [u0; v0];
end

function [pl,ql,pr,qr] = bcfun_neumann(xl,Ul,xr,Ur,t) %#ok<INUSD>
    pl = [0; 0];
    ql = [1; 1];
    pr = [0; 0];
    qr = [1; 1];
end

function [u0, v0] = initial_condition_gap(x,A,w,x_uedge,x_vedge)
    % smoothed Heaviside-type initial data
    H1 = 0.5 * (1 + tanh((x - x_uedge)/w));
    H2 = 0.5 * (1 + tanh((x - x_vedge)/w));

    u0 = A * (1 - H1);
    v0 = A * H2;
end

function [W0_at_0, W1_at_0] = leading_order_stefan_constants(a1)
% Constants in
% c ~ -gamma*W0(0)/(1 + gamma*W1(0)).

    U = linspace(0, 1, 10001);
    R = U .* (1-U) .* (U-a1);
    I = cumtrapz(U, R);
    tailI = max(I(end) - I, 0);
    W0 = -sqrt(2*tailI);

    W0_at_0 = W0(1);
    W1_at_0 = trapz(U, W0) / W0_at_0;
end

function [u0, v0] = initial_condition_stefan_small_branch(x,A,w,L0,a1,D)
% Contact initial data close to the corrected small-speed Stefan branch.
% The u-profile is the leading-order c=0 Stefan half-wave
% U_z = -sqrt(2*int_U^1 R(s;a1) ds), U(0)=0, U(-infty)=1.
% The v-profile is a smooth right-hand state in contact with this interface.

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

function Lx = refined_max_location(x, y, idx)
    if idx <= 1 || idx >= numel(x)
        Lx = x(idx);
        return;
    end

    x3 = x(idx-1:idx+1);
    y3 = y(idx-1:idx+1);

    if any(~isfinite(y3)) || y3(2) <= y3(1) || y3(2) <= y3(3)
        Lx = x(idx);
        return;
    end

    coeff = polyfit(x3, y3, 2);
    if coeff(1) >= 0
        Lx = x(idx);
        return;
    end

    Lx = -coeff(2) / (2*coeff(1));
    if Lx < x3(1) || Lx > x3(3)
        Lx = x(idx);
    end
end

function [bracket, scanTable] = find_finite_sign_bracket(D, bracket0, opts)
    nScan = pick(opts, 'n_bracket_scan', 15);
    targetC = pick(opts, 'target_c', 0);
    gScan = linspace(bracket0(1), bracket0(2), nScan);
    cScan = nan(size(gScan));

    for k = 1:numel(gScan)
        cScan(k) = run_one_case(D, gScan(k), opts).c;
    end

    finite = isfinite(cScan);
    scanTable = table(gScan(finite).', cScan(finite).', ...
        'VariableNames', {'gamma','c'});

    idxFinite = find(finite);
    bracket = [nan, nan];
    if numel(idxFinite) < 2
        return;
    end

    for k = 1:numel(idxFinite)-1
        i1 = idxFinite(k);
        i2 = idxFinite(k+1);
        r1 = cScan(i1) - targetC;
        r2 = cScan(i2) - targetC;
        if r1 == 0 || r1*r2 <= 0
            bracket = [gScan(i1), gScan(i2)];
            return;
        end
    end
end

function Lx = descending_level_location(x, urow, level)
    Lx = NaN;
    urow = urow(:);
    x = x(:);
    idx = find(urow(1:end-1) >= level & urow(2:end) < level, 1, 'last');
    if isempty(idx)
        return;
    end

    u1 = urow(idx);
    u2 = urow(idx+1);
    if ~isfinite(u1) || ~isfinite(u2) || abs(u2-u1) < eps
        Lx = x(idx);
        return;
    end

    theta = (level-u1)/(u2-u1);
    theta = max(0, min(1, theta));
    Lx = x(idx) + theta*(x(idx+1)-x(idx));
end

function Lx = find_interface_local_centroid(x, urow, vrow, prevL, ...
        overlapTol, halfWindow)
    q = urow .* vrow;
    q(~isfinite(q)) = 0;
    [qmax, idxMax] = max(q);

    if qmax <= overlapTol
        Lx = NaN;
        return;
    end

    if isfinite(prevL)
        nearPrevious = abs(x - prevL) <= 12;
        if nnz(nearPrevious) >= 3 && max(q(nearPrevious)) >= 0.05*qmax
            localIdx = find(nearPrevious);
            [~, k] = max(q(nearPrevious));
            idxMax = localIdx(k);
        end
    end

    i1 = max(1, idxMax - halfWindow);
    i2 = min(numel(x), idxMax + halfWindow);
    xLocal = x(i1:i2);
    qLocal = q(i1:i2);
    qIntegral = trapz(xLocal, qLocal);

    if ~isfinite(qIntegral) || qIntegral <= overlapTol
        Lx = NaN;
    else
        Lx = trapz(xLocal, xLocal .* qLocal) / qIntegral;
    end
end

function Lx = find_interface_reaction_peak(x, urow, vrow, prevL, active_tol)
% Track the reaction zone by the peak of u*v.
% This is more stable than u=v for very small D, where equality can be
% dominated by low-density tails after the populations have met.

    q = urow .* vrow;
    q(~isfinite(q)) = 0;

    if max(q) < active_tol^2
        Lx = NaN;
        return;
    end

    if ~isnan(prevL)
        windowWidth = 12;
        mask = abs(x - prevL) <= windowWidth;
        if nnz(mask) >= 3 && max(q(mask)) >= 0.05*max(q)
            idxLocal = find(mask);
            [~, loc] = max(q(mask));
            idx = idxLocal(loc);
        else
            [~, idx] = max(q);
        end
    else
        [~, idx] = max(q);
    end

    Lx = refined_max_location(x, q, idx);
end

function Lx = find_interface_uv_equal(x, urow, vrow, prevL, active_tol)
% Robust interface detection:
% 1) first look for sign changes of d = u-v
% 2) if none found, fall back to the point where |u-v| is minimal
%    near the previous interface location
% 3) reject only if both densities are too small

    d = urow - vrow;

    % ---------- Step 1: sign-change candidates ----------
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

        xi = x(i) + theta * (x(i+1) - x(i));
        ui = urow(i) + theta * (urow(i+1) - urow(i));
        vi = vrow(i) + theta * (vrow(i+1) - vrow(i));

        if min(ui,vi) < active_tol
            continue;
        end

        xcand(end+1,1) = xi; %#ok<AGROW>
        wcand(end+1,1) = ui + vi; %#ok<AGROW>
    end

    if ~isempty(xcand)
        if ~isnan(prevL)
            [~, jbest] = min(abs(xcand - prevL));
        else
            [~, jbest] = max(wcand);
        end
        Lx = xcand(jbest);
        return;
    end

    % ---------- Step 2: fallback to minimum |u-v| ----------
    absd = abs(d);

    if ~isnan(prevL)
        window = 20;
        [~, ic] = min(abs(x - prevL));
        i1 = max(1, ic - window);
        i2 = min(numel(x), ic + window);

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
