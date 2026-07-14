function plot_flux_matched_comparison()
% plot_flux_mu_asymptotics
%
% Compare numerical travelling-wave speeds with the small-c and near-critical
% asymptotic predictions for the reduced flux-type model.
%
% The comparison parameter is
%
%     eta = gamma/sqrt(D) = -1/(sqrt(3) W(0)).
%
% The small-c expansion gives
%
%     eta ~ -1/(sqrt(3) (W0(0) + c W1(0))).
%
% The near-critical expansion gives
%
%     eta ~ -1/(sqrt(3) psi0(0)) epsilon^(-p),
%     epsilon = c_* - c,    p = 1/(2a+1).

    clearvars -except ans;
    clc;
    close all;

    a = 0.30;
    cStar = sqrt(2) * (0.5 - a);
    F = @(U) U .* (1 - U) .* (U - a);

    cMin = -0.75 * cStar;
    cVals = [linspace(cMin, -1e-3, 120), ...
        linspace(1e-3, 0.98*cStar, 160), ...
        cStar - logspace(log10(0.02*cStar), -8, 100)];
    cVals = unique(cVals(:), 'stable');
    cVals = cVals(cVals < cStar);

    etaNum = nan(size(cVals));
    for i = 1:numel(cVals)
        Wfront = front_slope_from_phase_plane(F, cVals(i));
        etaNum(i) = -1 / (sqrt(3) * Wfront);
    end

    keep = isfinite(etaNum) & etaNum > 0 & isfinite(cVals);
    cVals = cVals(keep);
    etaNum = etaNum(keep);
    epsilonNum = cStar - cVals;

    [etaNum, idx] = sort(etaNum);
    cVals = cVals(idx);
    epsilonNum = epsilonNum(idx);

    [W00, W10] = small_c_endpoint_terms(a);
    etaSmallTargetMax = 10;
    cSmallMax = (-1 / (sqrt(3) * etaSmallTargetMax) - W00) / W10;
    cSmall = linspace(cMin, cSmallMax, 700).';
    etaSmall = -1 ./ (sqrt(3) * (W00 + cSmall .* W10));
    keepSmall = isfinite(etaSmall) & etaSmall > 0;
    cSmall = cSmall(keepSmall);
    etaSmall = etaSmall(keepSmall);

    p = 1 / (2*a + 1);
    sigma = -1 / sqrt(2);
    rho = 1 / (sqrt(2) * (2*a + 1));
    A = -gamma(2*a + 1) * gamma(3 - 2*a) / gamma(4);
    xi = (-sigma + rho)^rho * (-A)^(-sigma);
    psi0AtZero = -xi^(1 / (rho - sigma));

    epsilonGrid = logspace(log10(min(epsilonNum)), log10(max(epsilonNum)), 700).';
    cLarge = cStar - epsilonGrid;
    etaLarge = -1 ./ (sqrt(3) * psi0AtZero .* epsilonGrid.^p);
    keepLarge = isfinite(etaLarge) & etaLarge >= 5 & cLarge > 0;
    epsilonGrid = epsilonGrid(keepLarge);
    cLarge = cLarge(keepLarge);
    etaLarge = etaLarge(keepLarge);

    fitMask = epsilonNum > 1e-7 & epsilonNum < 1e-3;
    fitMaskRed = epsilonGrid > 1e-7 & epsilonGrid < 1e-3;
    fitBlue = polyfit(log10(epsilonNum(fitMask)), log10(etaNum(fitMask)), 1);
    fitRed = polyfit(log10(epsilonGrid(fitMaskRed)), log10(etaLarge(fitMaskRed)), 1);

    eta0 = -1 / (sqrt(3) * W00);
    etaJoin = 9;

    fprintf('Flux-type asymptotic comparison\n');
    fprintf('  a = %.6g\n', a);
    fprintf('  c_* = %.10f\n', cStar);
    fprintf('  W0(0) = %.10f\n', W00);
    fprintf('  W1(0) = %.10f\n', W10);
    fprintf('  zero-speed threshold eta0 = %.10f\n', eta0);
    fprintf('  p = %.10f\n', p);
    fprintf('  psi0(0) = %.10f\n', psi0AtZero);
    fprintf('  near-critical red slope = %.10f\n', fitRed(1));
    fprintf('  numerical blue slope = %.10f\n', fitBlue(1));
    fprintf('  eta range = [%.6e, %.6e]\n', min(etaNum), max(etaNum));
    fprintf('  c range = [%.6e, %.6e]\n', min(cVals), max(cVals));

    colNum = [0 0.2 1];
    colLarge = [1 0 0];
    colSmall = [0 0.7 0];

    fig1 = figure('Color', 'w', 'Position', [100 120 620 470]);

    ax1 = axes(fig1, 'Position', [0.12 0.14 0.80 0.78]);
    hold(ax1, 'on');
    smallPlotMask = etaSmall <= etaJoin;
    largePlotMask = etaLarge >= etaJoin;
    hNum1 = plot(ax1, etaNum, cVals, '-', 'Color', colNum, 'LineWidth', 1.8);
    hLarge1 = plot(ax1, etaLarge(largePlotMask), cLarge(largePlotMask), ...
        '--', 'Color', colLarge, 'LineWidth', 1.8);
    hSmall1 = plot(ax1, etaSmall(smallPlotMask), cSmall(smallPlotMask), ...
        '--', 'Color', colSmall, 'LineWidth', 1.6);
    yline(ax1, cStar, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
    annotation(fig1, 'textbox', [0.15 0.82 0.08 0.05], 'String', '$c^*$', ...
        'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 12);
    xlabel(ax1, '$\gamma/\sqrt{D}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(ax1, '$c$', 'Interpreter', 'latex', 'FontSize', 14);
    legend(ax1, [hNum1, hLarge1, hSmall1], ...
        {'Numerical solution', 'Large-$c$ expansion', 'Small-$c$ expansion'}, ...
        'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 8);
    xlim(ax1, [0, min(25, max(etaNum))]);
    ylim(ax1, [min(cVals), 1.04*cStar]);
    grid(ax1, 'on');
    box(ax1, 'on');
    set(ax1, 'FontSize', 11, 'TickLabelInterpreter', 'latex');

    fig2 = figure('Color', 'w', 'Position', [760 120 620 470]);

    ax2 = axes(fig2, 'Position', [0.12 0.14 0.80 0.78]);
    hold(ax2, 'on');
    posMask = cVals > 0 & etaNum > 0;
    posSmallMask = cSmall > 0 & etaSmall > 0 & etaSmall <= etaJoin;
    posLargeMask = cLarge > 0 & etaLarge > 0 & etaLarge >= etaJoin;
    hNum2 = loglog(ax2, etaNum(posMask), cVals(posMask), '-', ...
        'Color', colNum, 'LineWidth', 1.8);
    hLarge2 = loglog(ax2, etaLarge(posLargeMask), cLarge(posLargeMask), ...
        '--', 'Color', colLarge, 'LineWidth', 1.8);
    hSmall2 = loglog(ax2, etaSmall(posSmallMask), cSmall(posSmallMask), ...
        '--', 'Color', colSmall, 'LineWidth', 1.6);
    yline(ax2, cStar, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
    xlabel(ax2, '$\gamma/\sqrt{D}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(ax2, '$c$', 'Interpreter', 'latex', 'FontSize', 14);
    legend(ax2, [hNum2, hLarge2, hSmall2], ...
        {'Numerical solution', 'Large-$c$ expansion', 'Small-$c$ expansion'}, ...
        'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 9);
    grid(ax2, 'on');
    box(ax2, 'on');
    xlim(ax2, [min(etaNum(posMask)), max(etaNum(posMask))]);
    ylim(ax2, [min(cVals(posMask)), 1.04*cStar]);
    set(ax2, 'XScale', 'log', 'YScale', 'log');
    set(ax2, 'FontSize', 11, 'TickLabelInterpreter', 'latex');

    fig3 = figure('Color', 'w', 'Position', [420 650 620 470]);

    ax3 = axes(fig3, 'Position', [0.12 0.14 0.80 0.78]);
    hold(ax3, 'on');
    hNum3 = loglog(ax3, epsilonNum, etaNum, '-', 'Color', colNum, 'LineWidth', 1.8);
    hLarge3 = loglog(ax3, epsilonGrid, etaLarge, '--', ...
        'Color', colLarge, 'LineWidth', 1.8);
    xlabel(ax3, '$\varepsilon$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(ax3, '$\gamma/\sqrt{D}$', 'Interpreter', 'latex', 'FontSize', 14);
    legend(ax3, [hNum3, hLarge3], ...
        {'Numerical solution', 'Large-$c$ expansion'}, ...
        'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 9);
    grid(ax3, 'on');
    box(ax3, 'on');
    xlim(ax3, [min(epsilonNum), max(epsilonNum)]);
    ylim(ax3, [min(etaNum), max(etaNum)]);
    set(ax3, 'XScale', 'log', 'YScale', 'log');
    set(ax3, 'FontSize', 11, 'TickLabelInterpreter', 'latex');

end

function [W00, W10] = small_c_endpoint_terms(a)
    R = @(s) s .* (1 - s) .* (s - a);
    I0 = integral(R, 0, 1, 'RelTol', 1e-12, 'AbsTol', 1e-14);
    W00 = -sqrt(2 * I0);

    U = linspace(1e-7, 1 - 1e-7, 4000).';
    I = zeros(size(U));
    for i = 1:numel(U)
        I(i) = integral(R, U(i), 1, 'RelTol', 1e-11, 'AbsTol', 1e-13);
    end
    W0 = -sqrt(2 * max(I, 0));

    W0Integral = integral(@(s) interp1(U, W0, s, 'pchip', 'extrap'), ...
        0, 1, 'RelTol', 1e-10, 'AbsTol', 1e-12);
    W10 = W0Integral / W00;
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
