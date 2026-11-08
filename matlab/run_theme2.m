%RUN_THEME2 Stability of a discrete-time Takagi-Sugeno fuzzy system (Theme 2).
%   x(k+1) = h1(x) A1 x(k) + h2(x) A2 x(k),  h1 + h2 = 1,
%   A1 = [0.9 a; 0 0.8],  A2 = [0.9 0; a 0.8].
%   A. Values of a for which a common quadratic Lyapunov function
%      V = x'Px exists: P > 0, Ai'*P*Ai - P < 0 (i = 1, 2).
%   B. With h1(x) = min(max(x2 - 0.5, 0), 1): numerical study of the values
%      of a for which the origin is asymptotically stable, and an exact
%      condition for equilibria other than the origin.
%   Plain MATLAB. Writes ../results/theme2_results.txt and fig5_*, fig6_*.

here = fileparts(mfilename('fullpath'));
figDir = fullfile(here, '..', 'figures'); resDir = fullfile(here, '..', 'results');
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(resDir, 'dir'), mkdir(resDir); end
logFile = fullfile(resDir, 'theme2_results.txt');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);

mats = @(a) deal([0.9 a; 0 0.8], [0.9 0; a 0.8]);
D = diag([-1 1]);   % D*Ai(a)*D = Ai(-a) and x -> D*x keeps x2 (so h1): every result is symmetric in a

%% A. Common quadratic Lyapunov function
fprintf('A. Common quadratic Lyapunov function V = x''Px\n');
aGrid = 0:0.01:0.5;
fA = zeros(size(aGrid));
for k = 1:numel(aGrid)
    [A1, A2] = mats(aGrid(k));
    fA(k) = lyapunov_margin(A1, A2);
end
% Bisection for the boundary f(a) = 0.
lo = 0; hi = 0.5;
for it = 1:40
    mid = (lo + hi)/2;
    [A1, A2] = mats(mid);
    if lyapunov_margin(A1, A2) < 0, lo = mid; else, hi = mid; end
end
aLyap = lo;
fprintf('   A common P exists for |a| < %.4f (bisection; symmetric in a since D*Ai(a)*D = Ai(-a), D = diag(-1,1))\n', aLyap);
for a = [0.1 0.2 0.26]
    [A1, A2] = mats(a);
    [f, P] = lyapunov_margin(A1, A2);
    P = P/min(eig(P));                      % scale so that min eig(P) = 1
    fprintf('   a = %.2f: P = [%.4f %.4f; %.4f %.4f], max eig(Ai''PAi - P) = %.4f\n', a, P, ...
        max(max(eig(A1'*P*A1 - P)), max(eig(A2'*P*A2 - P))));
end
for a = [0.262 0.27 0.3]
    [A1, A2] = mats(a);
    fprintf('   a = %.3f: min over P of max eig = %+.2e (no common P)\n', a, lyapunov_margin(A1, A2));
end

%% B. The given membership function
fprintf('\nB. h1(x) = min(max(x2 - 0.5, 0), 1)\n');
aEq = sqrt(0.08);
fprintf('   Equilibria x ~= 0 need det(h1*A1 + h2*A2 - I) = 0.02 - h1*h2*a^2 = 0; h1*h2 <= 1/4,\n');
fprintf('   so they exist exactly when a^2 >= 0.08, i.e. |a| >= %.4f.\n', aEq);
fprintf('   Near the origin (x2 < 0.5) h1 = 0 and only A2 acts (eigenvalues 0.9, 0.8): the origin\n');
fprintf('   is locally asymptotically stable for every a.\n');

% Initial states: a 21 x 21 grid on [-10, 10]^2 and 1000 random states on [-100, 100]^2.
[g1, g2] = meshgrid(-10:1:10);
X0 = [g1(:)'; g2(:)'];
X0 = X0(:, any(X0 ~= 0, 1));
rng(1); X0 = [X0, 200*rand(2, 1000) - 100];
% Long-run outcome of each initial state after K steps:
%   1 origin (|x| < 1e-8), 2 equilibrium x ~= 0, 3 periodic orbit with
%   period 2..400, 4 bounded but not periodic up to 400 steps, 5 diverged.
K = 20000;
aSim = 0:0.0025:1;                  % results for a < 0 follow by symmetry
cls = zeros(5, numel(aSim));
for k = 1:numel(aSim)
    c = ts_classify(aSim(k), X0, K);
    cls(:, k) = histcounts(c, 0.5:1:5.5)';
end
frac = cls(1, :)/size(X0, 2);
fprintf('   Simulation: %d initial states (21 x 21 grid on [-10, 10]^2 and 1000 random on [-100, 100]^2),\n', size(X0, 2));
fprintf('   %d steps, a from 0 to 1 in steps of 0.0025. Outcome counts (origin / equilibrium / periodic /\n', K);
fprintf('   bounded non-periodic / diverged) where they change:\n');
prev = [];
for k = 1:numel(aSim)
    sig = cls(:, k) > 0;
    if ~isequal(sig, prev)
        fprintf('     from a = %.4f: %4d %4d %4d %4d %4d\n', aSim(k), cls(:, k));
        prev = sig;
    end
end
fprintf('   Every state reaches the origin for 0 <= a <= %.4f (and for the same |a| when a < 0).\n', ...
    aSim(find(frac < 1, 1) - 1));
fprintf('   No state diverged for any a in the scan (max over a of the diverged count: %d).\n', max(cls(5, :)));
for a = [0.282 0.2828 0.283 0.29 0.35 0.40 0.5 1]
    c = ts_classify(a, X0, K);
    fprintf('   a = %-6g: origin %4d, equilibrium %4d, periodic %4d, bounded non-periodic %4d, diverged %d', ...
        a, sum(c == 1), sum(c == 2), sum(c == 3), sum(c == 4), sum(c == 5));
    if a^2 > 0.08
        h = (1 + [-1 1]*sqrt(1 - 0.08/a^2))/2;   % roots of h1*(1 - h1) = 0.02/a^2
        fprintf('; predicted equilibria (%.3f, %.3f), (%.3f, %.3f)', [10*h.*a.*(0.5 + h); 0.5 + h]);
    end
    fprintf('\n');
end
[Xf, ~] = ts_iterate(0.29, X0, K);
h = (1 + [-1 1]*sqrt(1 - 0.08/0.29^2))/2; xe = [10*h.*0.29.*(0.5 + h); 0.5 + h];
st = vecnorm(Xf) >= 1e-8;
fprintf('   a = 0.29: states not at the origin end within %.1e of a predicted equilibrium\n', ...
    max(min([vecnorm(Xf(:, st) - xe(:, 1)); vecnorm(Xf(:, st) - xe(:, 2))], [], 1)));
diary('off');

%% Figures
f5 = figure('Color', 'w', 'Position', [100 100 1000 700]);
tl = tiledlayout(f5, 2, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
plot(ax, [-fliplr(aGrid) aGrid], [fliplr(fA) fA], 'LineWidth', 1.8);
yline(ax, 0, 'k-'); xline(ax, [-aLyap aLyap], '--', sprintf('|a| = %.4f', aLyap), 'LabelVerticalAlignment', 'bottom');
xlabel(ax, 'a'); ylabel(ax, 'min_P max_i \lambda_{max}(A_i^TPA_i - P)');
title(ax, 'A. Common quadratic Lyapunov function exists where < 0');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
aa = [-fliplr(aSim(2:end)) aSim];
plot(ax, aa, 100*[fliplr(frac(2:end)) frac], 'LineWidth', 1.8);
xline(ax, [-aLyap aLyap], '--', 'Color', [0.3 0.3 0.3]);
xline(ax, [-aEq aEq], '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
text(ax, aEq + 0.03, 60, '|a| = \surd0.08', 'Color', [0.85 0.33 0.1]);
xlabel(ax, 'a'); ylabel(ax, 'initial states reaching the origin (%)'); ylim(ax, [0 105]);
title(ax, 'B. Simulation with h_1 = min(max(x_2 - 0.5, 0), 1)');
% Outcome shares and a bifurcation diagram for a >= 0.
ax = nexttile(tl);
area(ax, aSim, 100*cls(1:4, :)'/size(X0, 2), 'LineStyle', 'none');
xlim(ax, [0.2 1]); ylim(ax, [0 100]); grid(ax, 'on');
xlabel(ax, 'a'); ylabel(ax, 'share of initial states (%)');
legend(ax, {'origin', 'equilibrium x \neq 0', 'periodic (2 to 400 steps)', 'bounded, not periodic'}, ...
    'Location', 'southeast', 'FontSize', 8);
title(ax, 'B. Long-run outcome (no state diverged)');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
aBif = 0.25:0.001:1;
for a = aBif
    [Xf, ~] = ts_iterate(a, [8 -8 4 9; 9 6 -9 3], 3000);
    Xf = Xf(:, vecnorm(Xf) > 1e-8);
    if isempty(Xf), plot(ax, a, 0, 'k.', 'MarkerSize', 3); continue; end
    xs = zeros(1, 150*size(Xf, 2));
    A1 = [0.9 a; 0 0.8]; A2 = [0.9 0; a 0.8];
    for n = 1:150
        h1 = min(max(Xf(2, :) - 0.5, 0), 1);
        Xf = h1.*(A1*Xf) + (1 - h1).*(A2*Xf);
        xs((n-1)*size(Xf, 2) + (1:size(Xf, 2))) = Xf(2, :);
    end
    plot(ax, a*ones(size(xs)), xs, 'k.', 'MarkerSize', 2);
end
xline(ax, aEq, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
xlabel(ax, 'a'); ylabel(ax, 'x_2 in the long run'); xlim(ax, [0.25 1]);
title(ax, 'B. Long-run x_2 for four starting states not attracted to 0');
exportgraphics(f5, fullfile(figDir, 'fig5_stability_vs_a.png'), 'Resolution', 150); close(f5);

f6 = figure('Color', 'w', 'Position', [100 100 1100 400]);
tl = tiledlayout(f6, 1, 3, 'TileSpacing', 'compact');
[s1, s2] = meshgrid(-10:4:10);
S0 = [s1(:)'; s2(:)']; S0 = S0(:, any(S0 ~= 0, 1));
for a = [0.25 0.35 0.5]
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    patch(ax, [-12 12 12 -12], [0.5 0.5 1.5 1.5], [1 0.85 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    for n = 1:size(S0, 2)
        Xt = ts_trajectory(a, S0(:, n), 300);
        plot(ax, Xt(1, :), Xt(2, :), '-', 'Color', [0 0.447 0.741 0.5], 'LineWidth', 0.8);
        plot(ax, Xt(1, 1), Xt(2, 1), 'o', 'MarkerSize', 3, 'Color', [0 0.447 0.741]);
    end
    plot(ax, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
    if a^2 > 0.08
        h = (1 + [-1 1]*sqrt(1 - 0.08/a^2))/2;
        plot(ax, 10*h.*a.*(0.5 + h), 0.5 + h, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    end
    axis(ax, [-12 12 -12 12]); xlabel(ax, 'x_1'); ylabel(ax, 'x_2');
    title(ax, sprintf('a = %g', a));
end
title(tl, 'Trajectories from a grid of initial states (shaded: 0.5 < x_2 < 1.5, both rules active; star: equilibrium)');
exportgraphics(f6, fullfile(figDir, 'fig6_phase_portraits.png'), 'Resolution', 150); close(f6);

function [X, Xprev] = ts_iterate(a, X, K)
% K steps of the fuzzy system for every column of X.
A1 = [0.9 a; 0 0.8]; A2 = [0.9 0; a 0.8];
Xprev = X;
for k = 1:K
    Xprev = X;
    h1 = min(max(X(2, :) - 0.5, 0), 1);
    X = h1.*(A1*X) + (1 - h1).*(A2*X);
end
end

function c = ts_classify(a, X0, K)
% Long-run outcome of each initial state (see the codes above).
[X, ~] = ts_iterate(a, X0, K);
c = zeros(1, size(X, 2));
c(vecnorm(X) < 1e-8) = 1;
c(c == 0 & ~all(isfinite(X), 1)) = 5;
c(c == 0 & vecnorm(X) > 1e6) = 5;
rest = find(c == 0);
Y = X(:, rest); Z = Y; per = zeros(1, numel(rest));
for p = 1:400
    [Z, ~] = ts_iterate(a, Z, 1);
    hit = per == 0 & vecnorm(Z - Y) < 1e-7*max(1, vecnorm(Y));
    per(hit) = p;
end
c(rest(per == 1)) = 2;
c(rest(per >= 2)) = 3;
c(rest(per == 0)) = 4;
end

function Xt = ts_trajectory(a, x0, K)
A1 = [0.9 a; 0 0.8]; A2 = [0.9 0; a 0.8];
Xt = zeros(2, K+1); Xt(:, 1) = x0;
for k = 1:K
    h1 = min(max(Xt(2, k) - 0.5, 0), 1);
    Xt(:, k+1) = h1*A1*Xt(:, k) + (1 - h1)*A2*Xt(:, k);
end
end
