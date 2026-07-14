function plot_stefan_corrected_kappa_asymptotics()
% plot_stefan_corrected_kappa_asymptotics
%
% Compare numerical travelling-wave speeds with the small-c and large-c
% asymptotic predictions for the reduced Stefan-type model.
%
% The small-c expansion uses gamma = kappa and
%     -c/kappa ~ W0(0) + c W1(0).
% The large-c expansion uses the inner-front scaling
%     kappa ~ -c/(psi0(0) epsilon^p), epsilon = c_* - c.

    clearvars -except ans;
    clc;
    close all;

    a = 0.30;
    cStar = sqrt(2) * (0.5 - a);
    F = @(U) U .* (1 - U) .* (U - a);

    cVals = [linspace(1e-3, 0.98*cStar, 160), ...
        cStar - logspace(log10(0.02*cStar), -8, 100)];
    cVals = unique(cVals(:), 'stable');
    cVals = cVals(cVals > 0 & cVals < cStar);

    kappaNum = nan(size(cVals));
    for i = 1:numel(cVals)
        Wfront = front_slope_from_phase_plane(F, cVals(i));
        kappaNum(i) = -cVals(i) / Wfront;
    end

    keep = isfinite(kappaNum) & kappaNum > 0 & isfinite(cVals);
    cVals = cVals(keep);
    kappaNum = kappaNum(keep);

    [kappaNum, idx] = sort(kappaNum);
    cVals = cVals(idx);
    epsilonNum = cStar - cVals;

    [W00, W10] = small_c_endpoint_terms(a);
    kappaSmall = linspace(0, 5, 700).';
    cSmall = -kappaSmall .* W00 ./ (1 + kappaSmall .* W10);
    keepSmall = isfinite(cSmall) & cSmall > 0 & cSmall < 1.2*cStar;
    kappaSmall = kappaSmall(keepSmall);
    cSmall = cSmall(keepSmall);

    p = 1 / (2*a + 1);
    sigma = -1 / sqrt(2);
    rho = 1 / (sqrt(2) * (2*a + 1));
    A = -gamma(2*a + 1) * gamma(3 - 2*a) / gamma(4);
    xi = (-sigma + rho)^rho * (-A)^(-sigma);
    psi0AtZero = -xi^(1 / (rho - sigma));

    kappaLarge = linspace(5, 25, 700).';
    epsilonLarge = (cStar ./ ((-psi0AtZero) .* kappaLarge)).^(1 / p);
    cLarge = cStar - epsilonLarge;
    cLarge(cLarge <= 0) = nan;

    epsilonGrid = logspace(log10(min(epsilonNum)), log10(max(epsilonNum)), 700).';
    cEps = cStar - epsilonGrid;
    kappaLargeEps = -cStar ./ (psi0AtZero .* epsilonGrid.^p);
    keepLargeEps = isfinite(kappaLargeEps) & kappaLargeEps > 0 & cEps > 0;
    epsilonGrid = epsilonGrid(keepLargeEps);
    kappaLargeEps = kappaLargeEps(keepLargeEps);

    fitMask = epsilonNum > 1e-7 & epsilonNum < 1e-3;
    fitMaskRed = epsilonGrid > 1e-7 & epsilonGrid < 1e-3;
    fitBlue = polyfit(log10(epsilonNum(fitMask)), log10(kappaNum(fitMask)), 1);
    fitRed = polyfit(log10(epsilonGrid(fitMaskRed)), log10(kappaLargeEps(fitMaskRed)), 1);

    fprintf('Stefan asymptotic comparison\n');
    fprintf('  a = %.6g\n', a);
    fprintf('  c_* = %.10f\n', cStar);
    fprintf('  W0(0) = %.10f\n', W00);
    fprintf('  W1(0) = %.10f\n', W10);
    fprintf('  p = %.10f\n', p);
    fprintf('  psi0(0) = %.10f\n', psi0AtZero);
    fprintf('  large-c red slope = %.10f\n', fitRed(1));
    fprintf('  numerical blue slope = %.10f\n', fitBlue(1));
    fprintf('  kappa range = [%.6e, %.6e]\n', min(kappaNum), max(kappaNum));

    colNum = [0 0.2 1];
    colLarge = [1 0 0];
    colSmall = [0 0.7 0];

    fig1 = figure('Color', 'w', 'Position', [100 120 620 470]);

    ax1 = axes(fig1, 'Position', [0.12 0.14 0.80 0.78]);
    hold(ax1, 'on');
    hNum1 = plot(ax1, kappaNum, cVals, '-', 'Color', colNum, 'LineWidth', 1.8);
    hLarge1 = plot(ax1, kappaLarge, cLarge, '--', 'Color', colLarge, 'LineWidth', 1.8);
    hSmall1 = plot(ax1, kappaSmall, cSmall, '--', 'Color', colSmall, 'LineWidth', 1.6);
    yline(ax1, cStar, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
    annotation(fig1, 'textbox', [0.15 0.82 0.08 0.05], 'String', '$c^*$', ...
        'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 12);
    xlabel(ax1, '$\kappa$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(ax1, '$c$', 'Interpreter', 'latex', 'FontSize', 14);
    legend(ax1, [hNum1, hLarge1, hSmall1], ...
        {'Numerical solution', 'Large-$c$ expansion', 'Small-$c$ expansion'}, ...
        'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 8);
    xlim(ax1, [0, min(25, max(kappaNum))]);
    ylim(ax1, [0, 1.04*cStar]);
    grid(ax1, 'on');
    box(ax1, 'on');
    set(ax1, 'FontSize', 11, 'TickLabelInterpreter', 'latex');

    fig2 = figure('Color', 'w', 'Position', [760 120 620 470]);

    ax2 = axes(fig2, 'Position', [0.12 0.14 0.80 0.78]);
    hold(ax2, 'on');
    posMask = kappaNum > 0 & cVals > 0;
    posSmallMask = kappaSmall > 0 & cSmall > 0;
    hNum2 = loglog(ax2, kappaNum(posMask), cVals(posMask), '-', ...
        'Color', colNum, 'LineWidth', 1.8);
    hLarge2 = loglog(ax2, kappaLarge, cLarge, '--', 'Color', colLarge, 'LineWidth', 1.8);
    hSmall2 = loglog(ax2, kappaSmall(posSmallMask), cSmall(posSmallMask), ...
        '--', 'Color', colSmall, 'LineWidth', 1.6);
    yline(ax2, cStar, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
    annotation(fig2, 'textbox', [0.15 0.82 0.08 0.05], 'String', '$c^*$', ...
        'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 12);
    xlabel(ax2, '$\kappa$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(ax2, '$c$', 'Interpreter', 'latex', 'FontSize', 14);
    legend(ax2, [hNum2, hLarge2, hSmall2], ...
        {'Numerical solution', 'Large-$c$ expansion', 'Small-$c$ expansion'}, ...
        'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 9);
    grid(ax2, 'on');
    box(ax2, 'on');
    xlim(ax2, [min(kappaNum(posMask)), max(kappaNum(posMask))]);
    ylim(ax2, [min(cVals(posMask)), 1.04*cStar]);
    set(ax2, 'XScale', 'log', 'YScale', 'log');
    set(ax2, 'FontSize', 11, 'TickLabelInterpreter', 'latex');

    fig3 = figure('Color', 'w', 'Position', [420 650 620 470]);

    ax3 = axes(fig3, 'Position', [0.12 0.14 0.80 0.78]);
    hold(ax3, 'on');
    epsPosMask = epsilonNum > 0 & kappaNum > 0;
    hNum3 = loglog(ax3, epsilonNum(epsPosMask), kappaNum(epsPosMask), ...
        '-', 'Color', colNum, 'LineWidth', 1.8);
    hLarge3 = loglog(ax3, epsilonGrid, kappaLargeEps, '--', ...
        'Color', colLarge, 'LineWidth', 1.8);
    xlabel(ax3, '$\varepsilon$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(ax3, '$\kappa$', 'Interpreter', 'latex', 'FontSize', 14);
    legend(ax3, [hNum3, hLarge3], ...
        {'Numerical solution', 'Large-$c$ expansion'}, ...
        'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 9);
    grid(ax3, 'on');
    box(ax3, 'on');
    xlim(ax3, [min(epsilonNum(epsPosMask)), max(epsilonNum(epsPosMask))]);
    ylim(ax3, [min(kappaNum(epsPosMask)), max(kappaNum(epsPosMask))]);
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
