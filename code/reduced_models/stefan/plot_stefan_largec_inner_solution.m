function plot_stefan_largec_inner_solution()
% plot_stefan_largec_inner_solution
%
% Inner solution for the leading-order Stefan travelling-wave problem in
% the large-c limit c -> c_*^-.

    clearvars -except ans;
    clc;
    close all;

    outDir = fileparts(mfilename('fullpath'));

    a = 0.30;
    p = 1 / (2*a + 1);

    sigma = -1 / sqrt(2);
    rho = 1 / (sqrt(2) * (2*a + 1));
    A = -gamma(2*a + 1) * gamma(3 - 2*a) / gamma(4);

    xi = (-sigma + rho)^rho * (-A)^(-sigma);
    psiAtZero = -xi^(1 / (rho - sigma));

    implicitFun = @(psi, phi) abs(psi - rho*phi).^rho .* ...
        abs(psi - sigma*phi).^(-sigma) - xi;

    phi = linspace(0, 2.25, 700).';
    psi = nan(size(phi));
    psi(1) = psiAtZero;

    for i = 2:numel(phi)
        psi(i) = continue_implicit_root(implicitFun, phi(i), psi(i-1), sigma);
    end

    psiFar1 = sigma * phi;
    psiFar2 = sigma * phi + A * phi.^(-p);
    psiFar2(phi == 0) = nan;

    fprintf('Large-c Stefan inner solution\n');
    fprintf('  a = %.6g\n', a);
    fprintf('  p = %.8f\n', p);
    fprintf('  sigma = %.8f\n', sigma);
    fprintf('  rho = %.8f\n', rho);
    fprintf('  A = %.8f\n', A);
    fprintf('  psi_0(0) = %.8f\n', psiAtZero);

    fig = figure('Color', 'w', 'Position', [120 120 720 480]);
    ax = axes(fig);
    hold(ax, 'on');

    plot(ax, phi, psi, '-', 'Color', [0 0 1], 'LineWidth', 2.0);
    plot(ax, phi, psiFar1, '--', 'Color', [0 0 0], 'LineWidth', 1.6);
    plot(ax, phi, psiFar2, '--', 'Color', [1 0 0], 'LineWidth', 1.8);

    xlabel(ax, '$\phi$', 'Interpreter', 'latex', 'FontSize', 17);
    ylabel(ax, '$\psi_0$', 'Interpreter', 'latex', 'FontSize', 17);
    legend(ax, {'$\psi_0(\phi)$', '$-\phi/\sqrt{2}$', ...
        '$-\phi/\sqrt{2}+A\phi^{-1/(2a+1)}$'}, ...
        'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 11);

    xlim(ax, [0 2.25]);
    ylim(ax, [-1.8 0]);
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', 13, 'LineWidth', 1.0, 'TickLabelInterpreter', 'latex');
end

function psi = continue_implicit_root(implicitFun, phi, psiPrev, sigma)
    localFun = @(psi) implicitFun(psi, phi);

    widths = [0.05 0.1 0.2 0.4 0.8 1.5];
    for width = widths
        bracket = [psiPrev - width, psiPrev + width];
        fBracket = localFun(bracket);

        if all(isfinite(fBracket)) && prod(fBracket) <= 0
            psi = fzero(localFun, bracket);
            return;
        end
    end

    guess = sigma * phi - 0.05;
    psi = fzero(localFun, guess);
end
