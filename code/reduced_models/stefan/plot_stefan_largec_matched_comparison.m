function plot_stefan_largec_matched_comparison()
% plot_stefan_largec_matched_comparison
%
% Matched large-c comparison for the leading-order Stefan travelling-wave
% problem. This follows the plotting strategy of the Allee-Stefan
% TWCompare.m script. This version keeps only the singular/singular
% comparison: the black curve is the outer expansion W0 + epsilon W1 and
% the blue front layer uses the singular inner matching branch.

    clearvars -except ans;
    clc;
    close all;

    outDir = fileparts(mfilename('fullpath'));

    a = 0.30;
    epsilon = 1e-3;
    cStar = sqrt(2) * (0.5 - a);

    F = @(U) U .* (1 - U) .* (U - a);

    p = 1 / (2*a + 1);
    sigma = -1 / sqrt(2);
    rho = 1 / (sqrt(2) * (2*a + 1));
    A = -gamma(2*a + 1) * gamma(3 - 2*a) / gamma(4);

    [zOuterRaw, UOuter, WOuter] = solve_outer_profile(F, cStar);
    c = cStar - epsilon;
    [UNum, WNum] = solve_combined_phase_plane(F, cStar, epsilon);

    zNum = z_from_phase_plane(UNum, WNum);
    zOuter = zOuterRaw;
    zOuter = align_at_U(UOuter, zOuter, UNum, zNum, 0.5);

    UFloorInner = 1e-5;
    UMatchProbe = 0.02;
    [UInnerProbe, WInnerProbe] = singular_inner_branch(a, epsilon, UMatchProbe, UFloorInner, sigma, A);
    UMatch = choose_matching_point(UNum, WNum, UInnerProbe, WInnerProbe, epsilon, p);
    [UInner, WInner] = singular_inner_branch(a, epsilon, UMatch, UFloorInner, sigma, A);

    zInner = z_from_phase_plane(UInner, WInner);
    zInner = align_at_U(UInner, zInner, UNum, zNum, UMatch);

    idxNum = UNum >= UMatch;
    idxInner = UInner < UMatch;

    UTwoOuter = UNum(idxNum);
    WTwoOuter = WNum(idxNum);
    zTwoOuter = zNum(idxNum);

    UTwoInner = UInner(idxInner);
    WTwoInner = WInner(idxInner);
    zTwoInner = zInner(idxInner);

    % TWAfterMatching-style profile alignment: use a common translation in z
    % so the numerical front is placed at z = 0 near the smallest plotted U.
    UShift = max([UFloorInner, min(UNum), min(UOuter), min(UInner)]);
    shiftCommon = -interp_monotone(UNum, zNum, UShift);

    fprintf('Large-c singular matched comparison\n');
    fprintf('  a = %.6g\n', a);
    fprintf('  epsilon = %.6g\n', epsilon);
    fprintf('  c_* = %.8f, c = %.8f\n', cStar, c);
    fprintf('  p = %.8f\n', p);
    fprintf('  U_match = %.8e\n', UMatch);
    fprintf('  U_floor_inner = %.8e\n', UFloorInner);
    fprintf('  z_shift = %.8e at U = %.8e\n', shiftCommon, UShift);
    fprintf('  W_num(U_match) = %.8e\n', interp_monotone(UNum, WNum, UMatch));
    fprintf('  W_inner(U_match) = %.8e\n', interp_monotone(UInner, WInner, UMatch));

    colTwo = [0 0.4470 0.7410];
    colNum = [0 0 0];
    colOne = [1 0 0];

    figProfile = figure('Color', 'w', 'Position', [120 120 720 480]);
    axProfile = axes(figProfile);
    hold(axProfile, 'on');

    hTwo = plot(axProfile, zTwoOuter + shiftCommon, UTwoOuter, '-', 'Color', colTwo, 'LineWidth', 2.0);
    plot(axProfile, zTwoInner + shiftCommon, UTwoInner, '-', 'Color', colTwo, 'LineWidth', 2.0);
    hNum = plot(axProfile, zNum + shiftCommon, UNum, '--', 'Color', colNum, 'LineWidth', 1.6);
    hOne = plot(axProfile, zOuter + shiftCommon, UOuter, '--', 'Color', colOne, 'LineWidth', 1.8);

    xlabel(axProfile, '$z$', 'Interpreter', 'latex', 'FontSize', 17);
    ylabel(axProfile, '$U$', 'Interpreter', 'latex', 'FontSize', 17);
    legend(axProfile, [hTwo, hNum, hOne], ...
        {'Two-term solution', 'Numerical solution', 'One-term solution'}, ...
        'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);
    xlim(axProfile, [-14 3]);
    ylim(axProfile, [-0.02 1.02]);
    grid(axProfile, 'on');
    box(axProfile, 'on');
    set(axProfile, 'FontSize', 13, 'LineWidth', 1.0, 'TickLabelInterpreter', 'latex');

    axInset = axes('Position', [0.26 0.28 0.25 0.25]);
    hold(axInset, 'on');
    plot(axInset, zTwoOuter + shiftCommon, UTwoOuter, '-', 'Color', colTwo, 'LineWidth', 1.4);
    plot(axInset, zTwoInner + shiftCommon, UTwoInner, '-', 'Color', colTwo, 'LineWidth', 1.4);
    plot(axInset, zNum + shiftCommon, UNum, '--', 'Color', colNum, 'LineWidth', 1.1);
    plot(axInset, zOuter + shiftCommon, UOuter, '--', 'Color', colOne, 'LineWidth', 1.2);
    xlim(axInset, [-16.2 0.2]);
    ylim(axInset, [0 max(0.03, 1.2 * UMatch)]);
    grid(axInset, 'on');
    box(axInset, 'on');
    set(axInset, 'FontSize', 8, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');

    figPhase = figure('Color', 'w', 'Position', [180 180 720 480]);
    axPhase = axes(figPhase);
    hold(axPhase, 'on');

    [UTwoOuterAsc, WTwoOuterAsc] = sort_by_U(UTwoOuter, WTwoOuter);
    [UTwoInnerAsc, WTwoInnerAsc] = sort_by_U(UTwoInner, WTwoInner);
    [UNumAsc, WNumAsc] = sort_by_U(UNum, WNum);
    [UOuterAsc, WOuterAsc] = sort_by_U(UOuter, WOuter);

    hTwoPhase = plot(axPhase, UTwoOuterAsc, WTwoOuterAsc, '-', 'Color', colTwo, 'LineWidth', 2.0);
    plot(axPhase, UTwoInnerAsc, WTwoInnerAsc, '-', 'Color', colTwo, 'LineWidth', 2.0);
    hNumPhase = plot(axPhase, UNumAsc, WNumAsc, '--', 'Color', colNum, 'LineWidth', 1.6);
    hOnePhase = plot(axPhase, UOuterAsc, WOuterAsc, '--', 'Color', colOne, 'LineWidth', 1.8);

    xlabel(axPhase, '$U$', 'Interpreter', 'latex', 'FontSize', 17);
    ylabel(axPhase, '$W$', 'Interpreter', 'latex', 'FontSize', 17);
    legend(axPhase, [hTwoPhase, hNumPhase, hOnePhase], ...
        {'Two-term solution', 'Numerical solution', 'One-term solution'}, ...
        'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);
    xlim(axPhase, [0 1]);
    grid(axPhase, 'on');
    box(axPhase, 'on');
    set(axPhase, 'FontSize', 13, 'LineWidth', 1.0, 'TickLabelInterpreter', 'latex');

    axPhaseInset = axes('Position', [0.26 0.24 0.28 0.25]);
    hold(axPhaseInset, 'on');
    plot(axPhaseInset, UTwoOuterAsc, WTwoOuterAsc, '-', 'Color', colTwo, 'LineWidth', 1.4);
    plot(axPhaseInset, UTwoInnerAsc, WTwoInnerAsc, '-', 'Color', colTwo, 'LineWidth', 1.4);
    plot(axPhaseInset, UNumAsc, WNumAsc, '--', 'Color', colNum, 'LineWidth', 1.1);
    plot(axPhaseInset, UOuterAsc, WOuterAsc, '--', 'Color', colOne, 'LineWidth', 1.2);
    xlim(axPhaseInset, [0, 1.2 * UMatch]);
    localMaskTwoOuter = UTwoOuterAsc <= 1.2 * UMatch;
    localMaskTwoInner = UTwoInnerAsc <= 1.2 * UMatch;
    localMaskNum = UNumAsc <= 1.2 * UMatch;
    localMaskOne = UOuterAsc <= 1.2 * UMatch;
    localVals = [WTwoOuterAsc(localMaskTwoOuter); WTwoInnerAsc(localMaskTwoInner); ...
        WNumAsc(localMaskNum); WOuterAsc(localMaskOne)];
    localVals = localVals(isfinite(localVals));
    if isempty(localVals)
        ylim(axPhaseInset, [-0.03 0.01]);
    else
        ylim(axPhaseInset, [1.15 * min(localVals), max(0.01, 1.15 * max(localVals))]);
    end
    grid(axPhaseInset, 'on');
    box(axPhaseInset, 'on');
    set(axPhaseInset, 'FontSize', 8, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');

end

function [z, U, W] = solve_outer_profile(F, c)
    rhs = @(~, Y) [Y(2); -c*Y(2) - F(Y(1))];
    zspan = linspace(-10, 50, 1600);
    Y0 = [1; -1e-3];
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'MaxStep', 1e-3);

    [z, Y] = ode15s(rhs, zspan, Y0, opts);
    U = Y(:,1);
    W = Y(:,2);

    keep = isfinite(U) & isfinite(W) & U >= -0.05 & U <= 1.05;
    z = z(keep);
    U = U(keep);
    W = W(keep);

    if U(1) < U(end)
        z = flipud(z);
        U = flipud(U);
        W = flipud(W);
    end
end

function [U, W] = solve_combined_phase_plane(F, cStar, epsilon)
    Uvals = linspace(1, 1e-8, 3000).';
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'MaxStep', 1e-3);

    W0start = -epsilon;

    dW0dU = @(U,W0) (-cStar.*W0 - F(U)) ./ W0;
    [U0, W0] = ode15s(dW0dU, Uvals, W0start, opts);

    U0 = U0(:);
    W0 = W0(:);
    W0Interp = @(Uq) interp_monotone(U0, W0, Uq);

    dW1dU = @(U,W1) 1 - ...
        (((-cStar.*W0Interp(U) - F(U)) ./ W0Interp(U).^2) + cStar ./ W0Interp(U)) .* W1;
    [U1, W1] = ode15s(dW1dU, Uvals, 0, opts);

    W1OnU0 = interp_monotone(U1, W1, U0);

    U = U0;
    W = W0 + epsilon * W1OnU0;

    keep = isfinite(U) & isfinite(W) & W < 0;
    U = U(keep);
    W = W(keep);
end

function [U, W] = singular_inner_branch(a, epsilon, Umatch, Ufloor, sigma, A)
    p = 1 / (2*a + 1);
    phiMatch = Umatch / epsilon^p;
    phiFloor = Ufloor / epsilon^p;
    phi = linspace(phiMatch, phiFloor, 3500).';
    psi = sigma .* phi + A .* phi.^(-p);

    U = epsilon^p * phi;
    W = epsilon^p * psi;

    keep = isfinite(U) & isfinite(W) & U > 0 & W < -1e-12;
    U = U(keep);
    W = W(keep);
end

function UMatch = choose_matching_point(UNum, WNum, UInner, WInner, epsilon, p)
    lo = max([0.5*epsilon, min(UInner), min(UNum)]);
    hi = min([max(UInner), max(UNum), 0.02]);
    candidates = logspace(log10(lo), log10(hi), 120).';

    mismatch = nan(size(candidates));
    for i = 1:numel(candidates)
        Wn = interp_monotone(UNum, WNum, candidates(i));
        Wi = interp_monotone(UInner, WInner, candidates(i));
        mismatch(i) = abs(Wn - Wi) / max(abs(Wn), 1e-12);
    end

    [~, idx] = min(mismatch);
    UMatch = candidates(idx);
end

function z = z_from_phase_plane(U, W)
    U = U(:);
    W = W(:);

    if U(1) > U(end)
        Udesc = U;
        Wdesc = W;
    else
        Udesc = flipud(U);
        Wdesc = flipud(W);
    end

    zdesc = cumtrapz(Udesc, 1 ./ Wdesc);
    zdesc = zdesc - interp1(Udesc, zdesc, min(Udesc), 'linear', 'extrap');

    z = interp1(Udesc, zdesc, U, 'linear', 'extrap');
end

function z = align_at_U(U, z, Uref, zref, Ualign)
    z = z + interp_monotone(Uref, zref, Ualign) - interp_monotone(U, z, Ualign);
end

function [Usort, Wsort] = sort_by_U(U, W)
    [Usort, idx] = sort(U(:), 'ascend');
    Wsort = W(idx);
end

function yq = interp_monotone(x, y, xq)
    x = x(:);
    y = y(:);
    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);

    [xSort, idx] = sort(x, 'ascend');
    ySort = y(idx);
    [xUnique, ia] = unique(xSort, 'stable');
    yUnique = ySort(ia);

    yq = interp1(xUnique, yUnique, xq, 'linear', 'extrap');
end

function psi = continue_implicit_root(implicitFun, phi, psiPrev, sigma)
    localFun = @(psi) implicitFun(psi, phi);
    widths = [0.02 0.05 0.1 0.2 0.4 0.8 1.6 3.2];

    for width = widths
        bracket = [psiPrev - width, psiPrev + width];
        fBracket = localFun(bracket);
        if all(isfinite(fBracket)) && prod(fBracket) <= 0
            psi = fzero(localFun, bracket);
            return;
        end
    end

    psi = fzero(localFun, sigma * phi - 0.05);
end
