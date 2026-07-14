function plot_stefan_smallc_asymptotics()
% plot_stefan_smallc_asymptotics
%
% Compare the numerical phase-plane trajectory for the leading-order
% Stefan travelling-wave problem with the small-c asymptotic expansion
%
%     W(U) = W0(U) + c W1(U) + ...
%
% where
%
%     W0(U) = -sqrt(2 int_U^1 R(s;a) ds),
%     W1(U) = int_U^1 W0(s) ds / W0(U).

    clearvars -except ans;
    clc;
    close all;

    outDir = fileparts(mfilename('fullpath'));

    a = 0.20;
    c = 0.02;

    Umin = 1e-5;
    Umax = 1 - 1e-6;
    nGrid = 2000;
    U = linspace(Umin, Umax, nGrid).';

    W_num = numerical_phase_plane(U, a, c);
    [W0, W1] = small_c_terms(U, a);
    W_two = W0 + c * W1;

    W_at_0 = endpoint_slope(a, c);
    gamma_stefan = -c / W_at_0;

    fprintf('Small-c Stefan phase-plane comparison\n');
    fprintf('  a = %.6g, c = %.6g\n', a, c);
    fprintf('  Numerical W(0) = %.8f\n', W_at_0);
    fprintf('  Stefan gamma = -c/W(0) = %.8f\n', gamma_stefan);

    z_num = travelling_wave_coordinate(U, W_num);
    z_one = travelling_wave_coordinate(U, W0);
    z_two = travelling_wave_coordinate(U, W_two);

    figProfile = figure('Color', 'w', 'Position', [120 120 720 480]);
    axProfile = axes(figProfile);
    hold(axProfile, 'on');

    plot(axProfile, z_num, U, '-', 'Color', [0 0 1], 'LineWidth', 2.0);
    plot(axProfile, z_one, U, '-', 'Color', [0 0 0], 'LineWidth', 2.0);
    plot(axProfile, z_two, U, '--', 'Color', [1 0 0], 'LineWidth', 2.0);

    xlabel(axProfile, '$z$', 'Interpreter', 'latex', 'FontSize', 17);
    ylabel(axProfile, '$U$', 'Interpreter', 'latex', 'FontSize', 17);
    legend(axProfile, {'Numerical (ODE)', 'One-term Asymptotic', 'Two-term Asymptotic'}, ...
        'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 11);

    grid(axProfile, 'on');
    box(axProfile, 'on');
    xlim(axProfile, [-8 0]);
    ylim(axProfile, [0 1]);
    set(axProfile, 'FontSize', 13, 'LineWidth', 1.0, 'TickLabelInterpreter', 'latex');

    figPhase = figure('Color', 'w', 'Position', [180 180 720 480]);
    axPhase = axes(figPhase);
    hold(axPhase, 'on');

    plot(axPhase, U, W_num, '-', 'Color', [0 0 1], 'LineWidth', 2.0);
    plot(axPhase, U, W0, '-', 'Color', [0 0 0], 'LineWidth', 2.0);
    plot(axPhase, U, W_two, '--', 'Color', [1 0 0], 'LineWidth', 2.0);

    xlabel(axPhase, '$U$', 'Interpreter', 'latex', 'FontSize', 17);
    ylabel(axPhase, '$W$', 'Interpreter', 'latex', 'FontSize', 17);
    legend(axPhase, {'Numerical (ODE)', 'One-term Asymptotic', 'Two-term Asymptotic'}, ...
        'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 11);

    grid(axPhase, 'on');
    box(axPhase, 'on');
    xlim(axPhase, [0 1]);
    set(axPhase, 'FontSize', 13, 'LineWidth', 1.0, 'TickLabelInterpreter', 'latex');
end

function z = travelling_wave_coordinate(U, W)
    z = cumtrapz(U, 1 ./ W);
    z = z - z(1);
end

function W = numerical_phase_plane(Uquery, a, c)
    F = @(U) U .* (1 - U) .* (U - a);

    Fp1 = a - 1;
    lambda = (-c + sqrt(c^2 - 4*Fp1)) / 2;
    U0 = 1 - 1e-7;
    W0 = lambda * (U0 - 1);

    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', 1e-3);
    [Usol, Wsol] = ode15s(@(U,W) (-c.*W - F(U)) ./ W, [U0, min(Uquery)], W0, opts);

    W = interp1(flipud(Usol), flipud(Wsol), Uquery, 'pchip', 'extrap');
end

function [W0, W1] = small_c_terms(U, a)
    R = @(s) s .* (1 - s) .* (s - a);

    I = zeros(size(U));
    for i = 1:numel(U)
        I(i) = integral(R, U(i), 1, 'RelTol', 1e-11, 'AbsTol', 1e-13);
    end

    W0 = -sqrt(2 * max(I, 0));

    J = zeros(size(U));
    for i = 1:numel(U)
        J(i) = integral(@(s) W0_interp(s, U, W0), U(i), 1, ...
            'RelTol', 1e-10, 'AbsTol', 1e-12);
    end

    W1 = J ./ W0;
end

function y = W0_interp(s, U, W0)
    y = interp1(U, W0, s, 'pchip', 'extrap');
end

function W0 = endpoint_slope(a, c)
    F = @(U) U .* (1 - U) .* (U - a);

    Fp1 = a - 1;
    lambda = (-c + sqrt(c^2 - 4*Fp1)) / 2;
    Ustart = 1 - 1e-7;
    Wstart = lambda * (Ustart - 1);

    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', 1e-3);
    [~, Wsol] = ode15s(@(U,W) (-c.*W - F(U)) ./ W, [Ustart, 1e-7], Wstart, opts);

    W0 = Wsol(end);
end
