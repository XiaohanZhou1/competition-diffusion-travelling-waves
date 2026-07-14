%% PDEPETwoinitial.m
% Full-domain coupled competition model with front tracking.
%
% PDE system on x in [0, Xmax]:
%   u_t = u_xx + R(u;a1) - (u v)/delta
%   v_t = D v_xx + S(v)  - gamma*(u v)/delta
%
% Initial idea:
%   - u starts on the left
%   - v starts on the right
%   - there is an initial gap between them
%
% We do NOT define L(t) initially.
% Instead:
%   1) track Xu(t): rightmost front of u
%   2) track Xv(t): leftmost front of v
%   3) define tc when Xu >= Xv
%   4) only for t >= tc, define L(t) by u(x,t)=v(x,t)
%
% Then estimate c from the linear phase of L(t) for t > 20.

clear; clc; close all;

%% ===== Parameters =====
D     = 10;        % try 2, 5, 10 ...
delta = 1e-3;      % start with 1e-2, then reduce later
gamma = 1;         % competition asymmetry
a1    = 0.1;       % Allee threshold for u
A     = 1.0;       % amplitude of initial Heaviside states

% reaction terms
R = @(u) u.*(1-u).*(u-a1);   % cubic strong Allee
S = @(v) v.*(1-v);           % logistic

%% ===== Domain and time =====
Xmax = 250;
Nx   = 1000;
x    = linspace(0, Xmax, Nx);

tEnd = 80;
Nt   = 320;
t    = linspace(0, tEnd, Nt);    % row vector for pdepe
tcol = t(:);                     % column vector for fitting etc.

m = 0;   % slab geometry for pdepe

%% ===== Initial conditions: separated smoothed Heaviside =====
gap = 80;
wIC = 1.0;

xmid    = 0.5*Xmax;
x_uedge = xmid - gap/2;   % u occupies x < x_uedge
x_vedge = xmid + gap/2;   % v occupies x > x_vedge

%% ===== Solve PDE =====
sol = pdepe(m, ...
    @(x,t,U,DUdx) pdefun(x,t,U,DUdx,D,delta,gamma,R,S), ...
    @(x) icfun_gap(x,A,wIC,x_uedge,x_vedge), ...
    @bcfun_neumann, ...
    x, t);

u = sol(:,:,1);   % size Nt x Nx
v = sol(:,:,2);   % size Nt x Nx

%% ===== Track fronts Xu(t), Xv(t) BEFORE interface forms =====
% Thresholds for fronts
eta_u = 0.01;
eta_v = 0.01;

Xu = nan(Nt,1);   % rightmost x such that u >= eta_u
Xv = nan(Nt,1);   % leftmost  x such that v >= eta_v

for k = 1:Nt
    idx_u = find(u(k,:) >= eta_u, 1, 'last');
    if ~isempty(idx_u)
        Xu(k) = x(idx_u);
    end

    idx_v = find(v(k,:) >= eta_v, 1, 'first');
    if ~isempty(idx_v)
        Xv(k) = x(idx_v);
    end
end

%% ===== Detect contact time tc =====
tc_idx = find(Xu >= Xv, 1, 'first');

if isempty(tc_idx)
    tc = NaN;
    fprintf('No contact detected within the simulated time interval.\n');
else
    tc = t(tc_idx);
    fprintf('Contact detected at approximately tc = %.6f\n', tc);
end

%% ===== Define interface L(t) ONLY AFTER contact =====
% New definition:
%   L(t) = x where u(x,t) = v(x,t)
%
% We detect sign changes of d(x)=u-v and linearly interpolate.
% To avoid picking trivial crossings where u≈v≈0, require both to exceed active_tol.

L = nan(Nt,1);
active_tol = 1e-5;

if ~isempty(tc_idx)
    prevL = NaN;

    for k = tc_idx:Nt
        Lk = find_interface_uv_equal(x, u(k,:), v(k,:), prevL, active_tol);
        L(k) = Lk;

        if ~isnan(Lk)
            prevL = Lk;
        end
    end
end

%% ===== Estimate wave speed c from linear phase of L(t), for t > 20 =====
t_fit_start = 20;
fit_mask = (tcol > t_fit_start) & ~isnan(L);

if nnz(fit_mask) >= 2
    p = polyfit(tcol(fit_mask), L(fit_mask), 1);
    slopeL = p(1);
    interceptL = p(2);

    % If L(t) decreases linearly, slope is negative.
    % The propagation speed c is taken as the positive magnitude.
    c_est = abs(slopeL);

    L_fit = polyval(p, tcol(fit_mask));

    % R^2
    ydata = L(fit_mask);
    SSres = sum((ydata - L_fit).^2);
    SStot = sum((ydata - mean(ydata)).^2);

    if SStot > 0
        R2 = 1 - SSres/SStot;
    else
        R2 = NaN;
    end

    fprintf('Linear fit of L(t) for t > %.2f:\n', t_fit_start);
    fprintf('    L(t) ~ %.8f * t + %.8f\n', slopeL, interceptL);
    fprintf('Estimated speed c = |slope| = %.8f\n', c_est);
    fprintf('R^2 = %.8f\n', R2);
else
    slopeL     = NaN;
    interceptL = NaN;
    c_est      = NaN;
    R2         = NaN;
    fprintf('Not enough valid L(t) points for linear fitting when t > %.2f.\n', t_fit_start);
end

%% ===== Plot 1: heatmap of u =====
figure;
imagesc(x, t, u);
axis xy;
colorbar;
xlabel('x');
ylabel('t');
title('u(x,t)');

%% ===== Plot 2: heatmap of v =====
figure;
imagesc(x, t, v);
axis xy;
colorbar;
xlabel('x');
ylabel('t');
title('v(x,t)');

%% ===== Plot 3: front locations Xu(t), Xv(t) =====
figure; hold on;
plot(t, Xu, 'LineWidth', 1.5, 'DisplayName', 'X_u(t): right front of u');
plot(t, Xv, 'LineWidth', 1.5, 'DisplayName', 'X_v(t): left front of v');
xlabel('t');
ylabel('front location');
box on;
grid on;
legend('Location','best');
title('Tracked fronts');

%% ===== Plot 4: interface L(t) after contact =====
figure; hold on;

plot(t, Xu, '--', 'LineWidth', 1.2, 'DisplayName', 'X_u(t)');
plot(t, Xv, '--', 'LineWidth', 1.2, 'DisplayName', 'X_v(t)');

if ~all(isnan(L))
    plot(t, L, 'k', 'LineWidth', 1.8, 'DisplayName', 'L(t): u=v');
end

if nnz(fit_mask) >= 2
    tt_fit = tcol(fit_mask);
    plot(tt_fit, polyval(p, tt_fit), 'r-', 'LineWidth', 2.0, ...
        'DisplayName', sprintf('fit: slope=%.4f, c=%.4f', slopeL, c_est));
end

xlabel('t');
ylabel('location');
box on;
grid on;
legend('Location','best');
title('Interface tracking using u=v');

%% ===== Plot 5: u and v profiles on the same figure =====
figure; hold on;

tshow = [0 10 20 30 40 50 60];
for j = 1:numel(tshow)
    [~,k] = min(abs(t - tshow(j)));

    plot(x, u(k,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('u, t=%.0f', t(k)));
    plot(x, v(k,:), '--', 'LineWidth', 1.3, ...
        'DisplayName', sprintf('v, t=%.0f', t(k)));

    if ~isnan(L(k))
        xline(L(k), ':', 'HandleVisibility', 'off');
    end
end

xlabel('x');
ylabel('density');
grid on;
ylim([-0.05, 1.05]);
legend('Location','best');
title('Profiles and interface position');

disp('Done.');

%% ===== Local functions =====
function [c,f,s] = pdefun(x,t,U,DUdx,D,delta,gamma,R,S) %#ok<INUSD>
    u = U(1);
    v = U(2);

    ux = DUdx(1);
    vx = DUdx(2);

    c = [1; 1];
    f = [ux; D*vx];
    s = [R(u) - (u*v)/delta;
         S(v) - gamma*(u*v)/delta];
end

function U0 = icfun_gap(x,A,w,x_uedge,x_vedge)
    % Smooth Heaviside approximation
    % H(z) ~ 0.5*(1+tanh(z/w))
    %
    % u ~ A for x < x_uedge, 0 otherwise
    % v ~ A for x > x_vedge, 0 otherwise

    H1 = 0.5*(1 + tanh((x - x_uedge)/w));
    H2 = 0.5*(1 + tanh((x - x_vedge)/w));

    u0 = A*(1 - H1);
    v0 = A*H2;

    U0 = [u0; v0];
end

function [pl,ql,pr,qr] = bcfun_neumann(xl,Ul,xr,Ur,t) %#ok<INUSD>
    % Neumann BC at both ends:
    % u_x = 0, v_x = 0
    %
    % p + q*f = 0
    % Here f = [u_x; D*v_x], so q=1 and p=0 imply
    % u_x = 0 and D*v_x = 0 => v_x = 0

    pl = [0; 0];
    ql = [1; 1];
    pr = [0; 0];
    qr = [1; 1];
end

function Lx = find_interface_uv_equal(x, urow, vrow, prevL, active_tol)
    % Find interface location Lx where u=v by locating zeros of d=u-v
    % and using linear interpolation.
    %
    % If multiple candidates exist:
    %   - choose the one nearest prevL if prevL is known
    %   - otherwise choose the one with largest overlap weight
    %
    % Crossings where both u and v are too small are rejected.

    d = urow - vrow;

    % sign change candidates
    idx = find(d(1:end-1).*d(2:end) <= 0);

    if isempty(idx)
        Lx = NaN;
        return;
    end

    xcand = [];
    wcand = [];

    for j = 1:numel(idx)
        i = idx(j);

        d1 = d(i);
        d2 = d(i+1);

        % linear interpolation
        if abs(d2 - d1) < 1e-14
            theta = 0.5;
        else
            theta = -d1/(d2 - d1);
        end

        if theta < 0 || theta > 1
            continue;
        end

        xi = x(i) + theta*(x(i+1) - x(i));
        ui = urow(i) + theta*(urow(i+1) - urow(i));
        vi = vrow(i) + theta*(vrow(i+1) - vrow(i));

        % reject trivial crossing in near-empty region
        if min(ui, vi) < active_tol
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
        [~,jbest] = min(abs(xcand - prevL));
    else
        [~,jbest] = max(wcand);
    end

    Lx = xcand(jbest);
end
