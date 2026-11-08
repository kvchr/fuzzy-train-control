%RUN_THEME1 Fuzzy speed and stopping control of an electric train (Theme 1).
%   A. Nominal trip with the default controller.
%   B. How the closed loop depends on the controller: membership functions
%      (number and shape of the sets), rule consequents and inference
%      (defuzzification).
%   C. Robustness: mass +-10 % and a track slope of +-2 %.
%   Needs the Fuzzy Logic Toolbox. Writes ../results/theme1_results.txt and
%   the figures fig1_* to fig4_* in ../figures.

here = fileparts(mfilename('fullpath'));
figDir = fullfile(here, '..', 'figures'); resDir = fullfile(here, '..', 'results');
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(resDir, 'dir'), mkdir(resDir); end
logFile = fullfile(resDir, 'theme1_results.txt');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);

P = train_params();
[fis, C] = train_fis();
fprintf('Theme 1: train m = %g kg, a1 = %g, a2 = %g, a3 = %g, b = %g; trip L = %g km, Vref = %g m/s (%g km/h)\n', ...
    P.m, P.a1, P.a2, P.a3, P.b, P.L/1000, P.Vref, 3.6*P.Vref);
fprintf(['Controller: %d rules (%s); braking-zone centres d_i = [%s] m for speeds v_i = [%s] m/s\n' ...
    '            (braking profile adec = %g m/s^2); K = %g, Kv = %g, uff(Vref) = %.3f; triangles, wtaver\n'], ...
    numel(C.names), strjoin(C.names, ', '), num2str(C.dNodes), num2str(C.vNodes), C.adec, C.K, C.Kv, C.uffRef);
fprintf('Near the station: u = Kp*d - K*v with Kp = %.3f 1/m; damping ratio zeta = %.3f (90 kg), %.3f (100 kg), %.3f (110 kg)\n', ...
    C.Kp, C.zeta(90), C.zeta(100), C.zeta(110));
[Dg, Vg] = meshgrid(0:1:400, 0:0.05:22);
Ug = zeros(size(Dg));
for k = 1:numel(Dg), Ug(k) = min(max(evalfis(fis, [Dg(k), Vg(k)]), -1), 1); end
[gd, gv] = gradient(Ug, 1, 0.05);
fprintf('Steepness of u(d, v) for 0 <= d <= 400 m, 0 <= v <= 22 m/s: max |du/dv| = %.2f s/m, max |du/dd| = %.3f 1/m\n', ...
    max(abs(gv(:))), max(abs(gd(:))));

% Reference: fastest trip allowed by |u| <= 1 and v <= Vref (full power,
% cruise, full brake), to judge what "fast" means for this train.
fprintf('Fastest trip with |u| <= 1 and v <= Vref (full power, cruise, full brake): %.1f s\n\n', fastest_trip(P));

fprintf(['Columns: tArr arrival time (|L - x| <= 0.5 m and |v| <= 0.05 m/s from then on; - = never), ' ...
    'over = max(x) - L if positive,\n         stop = final x - L, vCr mean speed for 5 km < x < 25 km, ' ...
    'vMax max speed, aMin/aMax extreme acceleration (m/s^2),\n         |du| largest change of u in one ' ...
    '0.1 s sample, zeta damping ratio of the linearized stop (m = 100 kg)\n\n']);
z = @(x) x.*(abs(x) >= 5e-4);               % print 0 instead of -0
tstr = @(t) char(string(sprintf('%8.1f', t)).replace('     NaN', '       -'));
hdr = sprintf('%-32s %8s %8s %8s %7s %7s %6s %6s %6s %6s', 'case', 'tArr(s)', 'over(m)', 'stop(m)', ...
    'vCr', 'vMax', 'aMin', 'aMax', '|du|', 'zeta');
row = @(name, R, zt) fprintf('%-32s %s %8.3f %8.3f %7.3f %7.3f %6.2f %6.2f %6.3f %6.2f\n', ...
    name, tstr(R.tArrive), z(max(R.overshoot, 0)), z(R.stopErr), R.vCruise, R.vMax, R.aMin, R.aMax, R.duMax, zt);

%% A. Nominal trip
R0 = simulate_train(fis, P);
fprintf('A. Nominal trip\n%s\n', hdr); row('default controller', R0, C.zeta(100));
fprintf('   speed never negative: %d, position never decreases: %d, min distance to go %.4f m\n\n', ...
    all(R0.v >= -1e-9), all(diff(R0.x) >= -1e-9), min(P.L - R0.x));

%% B. Dependence on the controller settings
fprintf('B. Dependence on the controller (nominal train)\n%s\n', hdr);
cases = {
    'default controller',             {}
    'MFs: 2 braking-zone steps',      {'nNodes', 2}
    'MFs: 10 braking-zone steps',     {'nNodes', 10}
    'MFs: Gaussian, wtaver',          {'shape', 'gauss'}
    'MFs: Gaussian, wtsum',           {'shape', 'gauss', 'defuzz', 'wtsum'}
    'MFs: triangles, wtsum',          {'defuzz', 'wtsum'}
    'rules: K = 0.5',                 {'K', 0.5}
    'rules: K = 2',                   {'K', 2}
    'rules: adec = 0.5 m/s^2',        {'adec', 0.5}
    'rules: adec = 2 m/s^2',          {'adec', 2}
    'rules: no feed-forward',         {'ff', false}
    'rules: Kv = 0.25',               {'Kv', 0.25}
    'rules: Kv = 0.05',               {'Kv', 0.05}
    };
RB = cell(size(cases, 1), 1);
for i = 1:size(cases, 1)
    [fi, Ci] = train_fis(cases{i, 2}{:});
    RB{i} = simulate_train(fi, P);
    row(cases{i, 1}, RB{i}, Ci.zeta(100));
end
fg = train_fis('shape', 'gauss'); dd = linspace(-20, 400, 4201);
sg = zeros(size(dd));
for k = 1:numel(fg.Inputs(1).MembershipFunctions), sg = sg + evalmf(fg.Inputs(1).MembershipFunctions(k), dd); end
fprintf('Gaussian MFs: the degrees add up to between %.3f and %.3f (triangles: exactly 1)\n', min(sg), max(sg));
fprintf('\nB2. Cruise rule with the heavier train (m = 110 kg, level track)\n%s\n', hdr);
P110 = train_params('m', 110);
for c = {{'Kv = 1.25 (default)', {}}, {'no feed-forward', {'ff', false}}, {'Kv = 0.25', {'Kv', 0.25}}, {'Kv = 0.05', {'Kv', 0.05}}}
    [fi, Ci] = train_fis(c{1}{2}{:});
    row(c{1}{1}, simulate_train(fi, P110), Ci.zeta(110));
end
fprintf('\nB3. Braking-zone gain on a 2 %% downhill grade (m = 110 kg)\n%s\n', hdr);
Pdn = train_params('m', 110, 'grade', -0.02);
for c = {{'K = 0.5', {'K', 0.5}}, {'K = 1.25 (default)', {}}, {'K = 2', {'K', 2}}}
    [fi, Ci] = train_fis(c{1}{2}{:});
    row(c{1}{1}, simulate_train(fi, Pdn), Ci.zeta(110));
end
fprintf('\n');

%% C. Robustness: mass +-10 % and slope
fprintf('C. Robustness (default controller); grade > 0 is uphill\n%s\n', hdr);
masses = [90 100 110]; grades = [-0.02 0 0.02];
RC = cell(numel(masses), numel(grades));
for i = 1:numel(masses)
    for j = 1:numel(grades)
        RC{i, j} = simulate_train(fis, train_params('m', masses(i), 'grade', grades(j)));
        row(sprintf('m = %d kg, grade %+g %%', masses(i), 100*grades(j)), RC{i, j}, C.zeta(masses(i)));
    end
end
gMax = (P.b - (P.a1 + 110*P.a2)*P.Vref - P.a3*P.Vref^2)/(110*P.g);
fprintf('Steepest uphill grade at which u = 1 still holds Vref with m = 110 kg: %.2f %%\n', 100*tan(asin(gMax)));
fprintf('Predicted stopping offset m*g*sin(theta)/(b*Kp): %.3f m (m = 110 kg, 2 %% grade)\n', ...
    110*P.g*sin(atan(0.02))/(P.b*C.Kp));
diary('off');

%% Figures
cols = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556; 0.3 0.3 0.3];

% Fig 1: membership functions and control surface.
f1 = figure('Color', 'w', 'Position', [100 100 1050 400]);
tl = tiledlayout(f1, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
dd = linspace(-20, 260, 2801);
nm = numel(fis.Inputs(1).MembershipFunctions);
cm = lines(nm);
for k = 1:nm
    plot(ax, dd, evalmf(fis.Inputs(1).MembershipFunctions(k), dd), 'LineWidth', 1.8, 'Color', cm(k, :));
end
xlabel(ax, 'distance to the station d (m)'); ylabel(ax, 'membership degree');
legend(ax, C.names, 'Location', 'eastoutside'); ylim(ax, [-0.05 1.1]); xlim(ax, [-20 260]);
title(ax, 'Fuzzy sets of d');
ax = nexttile(tl); hold(ax, 'on');
w = Dg(1, :) <= 260;
contourf(ax, Dg(:, w), Vg(:, w), Ug(:, w), -1:0.1:1, 'LineStyle', 'none');
colormap(ax, turbo); clim(ax, [-1 1]); cb = colorbar(ax); cb.Label.String = 'u (clipped to [-1, 1])';
dp = linspace(0, 260, 1041);
vHat = zeros(size(dp));
for k = 1:numel(dp)
    mu = arrayfun(@(mf) evalmf(mf, dp(k)), fis.Inputs(1).MembershipFunctions);
    vHat(k) = sum(mu(2:end).*C.vNodes)/sum(mu(2:end));
end
h1 = plot(ax, dp, min(C.Vref, sqrt(2*C.adec*dp)), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.6);
h2 = plot(ax, dp, vHat, 'k-', 'LineWidth', 1.4);
h3 = plot(ax, C.dNodes, C.vNodes, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
legend(ax, [h1 h2 h3], {'min(V_{ref}, (2a_{dec}d)^{1/2})', 'target speed of the rules', 'centres (d_i, v_i)'}, ...
    'Location', 'southeast', 'FontSize', 8);
xlabel(ax, 'distance to the station d (m)'); ylabel(ax, 'speed v (m/s)'); xlim(ax, [0 260]); ylim(ax, [0 22]);
title(ax, 'Control surface u(d, v) near the station');
exportgraphics(f1, fullfile(figDir, 'fig1_controller.png'), 'Resolution', 150); close(f1);

% Fig 2: nominal trip.
f2 = figure('Color', 'w', 'Position', [100 100 1000 620]);
tl = tiledlayout(f2, 3, 2, 'TileSpacing', 'compact', 'TileIndexing', 'columnmajor');
tA = R0.t(find(P.L - R0.x < 300, 1)) - 5;   % start of the arrival window
ax = nexttile(tl); plot(ax, R0.t, R0.x/1000, 'LineWidth', 1.6); grid(ax, 'on'); ylabel(ax, 'x (km)');
title(ax, 'Whole trip'); xline(ax, tA, ':');
ax = nexttile(tl); plot(ax, R0.t, R0.v, 'LineWidth', 1.6); grid(ax, 'on'); ylabel(ax, 'v (m/s)'); ylim(ax, [0 22]);
yline(ax, P.Vref, '--', 'V_{ref}', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom'); xline(ax, tA, ':');
ax = nexttile(tl); plot(ax, R0.t, R0.u, 'LineWidth', 1.6); grid(ax, 'on'); ylabel(ax, 'u'); ylim(ax, [-1.1 1.1]);
xlabel(ax, 'time t (s)'); xline(ax, tA, ':');
w = R0.t >= tA & R0.t <= R0.tArrive + 15;
ax = nexttile(tl); plot(ax, R0.t(w), P.L - R0.x(w), 'LineWidth', 1.6); grid(ax, 'on'); ylabel(ax, 'L - x (m)');
title(ax, 'Arrival (zoom)'); xline(ax, R0.tArrive, '--', 'arrival');
ax = nexttile(tl); plot(ax, R0.t(w), R0.v(w), 'LineWidth', 1.6); grid(ax, 'on'); ylabel(ax, 'v (m/s)');
ax = nexttile(tl); plot(ax, R0.t(w), R0.u(w), 'LineWidth', 1.6); grid(ax, 'on'); ylabel(ax, 'u'); xlabel(ax, 'time t (s)');
exportgraphics(f2, fullfile(figDir, 'fig2_nominal_trip.png'), 'Resolution', 150); close(f2);

% Fig 3: parameter study.
f3 = figure('Color', 'w', 'Position', [100 100 1000 640]);
tl = tiledlayout(f3, 2, 2, 'TileSpacing', 'compact');
groups = {
    'Fuzzy sets: number and shape',            [1 2 3 4],     'dist'
    'Braking-zone rules: K and a_{dec}',       [1 7 8 9 10],  'dist'
    'Last metres (overshoot when \zeta < 1)',  [1 3 4 7 10],  'zoom'
    'Cruise rule: K_v and feed-forward',       [1 11 12 13],  'time'
    };
for g = 1:size(groups, 1)
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    idx = groups{g, 2};
    for n = 1:numel(idx)
        Rn = RB{idx(n)};
        if strcmp(groups{g, 3}, 'time')
            plot(ax, Rn.t, Rn.v, 'LineWidth', 1.5, 'Color', cols(min(n, 5), :));
        else
            plot(ax, P.L - Rn.x, Rn.v, 'LineWidth', 1.5, 'Color', cols(min(n, 5), :));
        end
    end
    switch groups{g, 3}
        case 'dist'
            set(ax, 'XDir', 'reverse'); xlim(ax, [0 300]); ylim(ax, [0 22]); xlabel(ax, 'distance to go L - x (m)');
        case 'zoom'
            set(ax, 'XDir', 'reverse'); xlim(ax, [-1.5 8]); ylim(ax, [0 4]); xline(ax, 0, 'k-');
            xlabel(ax, 'distance to go L - x (m)');
        case 'time'
            xlim(ax, [0 120]); ylim(ax, [0 22]); xlabel(ax, 'time t (s)');
    end
    ylabel(ax, 'v (m/s)'); title(ax, groups{g, 1});
    legend(ax, cases(idx, 1), 'Location', 'southwest', 'FontSize', 7);
end
exportgraphics(f3, fullfile(figDir, 'fig3_parameter_study.png'), 'Resolution', 150); close(f3);

% Fig 4: robustness.
f4 = figure('Color', 'w', 'Position', [100 100 1000 380]);
tl = tiledlayout(f4, 1, 2, 'TileSpacing', 'compact');
ax1 = nexttile(tl); hold(ax1, 'on'); grid(ax1, 'on');
ax2 = nexttile(tl); hold(ax2, 'on'); grid(ax2, 'on');
ls = {'--', '-', ':'}; names = {};
for i = 1:numel(masses)
    for j = 1:numel(grades)
        Rij = RC{i, j};
        plot(ax1, Rij.t, Rij.v, ls{j}, 'LineWidth', 1.3, 'Color', cols(i, :));
        plot(ax2, P.L - Rij.x, Rij.v, ls{j}, 'LineWidth', 1.3, 'Color', cols(i, :));
        names{end+1} = sprintf('%d kg, %+g %%', masses(i), 100*grades(j)); %#ok<SAGROW>
    end
end
xlim(ax1, [0 80]); ylim(ax1, [0 22]); xlabel(ax1, 'time t (s)'); ylabel(ax1, 'v (m/s)'); title(ax1, 'Start: speed');
set(ax2, 'XDir', 'reverse'); xlim(ax2, [-0.5 4]); ylim(ax2, [0 2.5]);
xline(ax2, 0, 'k-'); xlabel(ax2, 'distance to go L - x (m)'); ylabel(ax2, 'v (m/s)'); title(ax2, 'Last metres before the stop');
legend(ax2, names, 'Location', 'eastoutside', 'FontSize', 8);
exportgraphics(f4, fullfile(figDir, 'fig4_robustness.png'), 'Resolution', 150); close(f4);

function tMin = fastest_trip(P)
% Full power until Vref, cruise at Vref, full brake to rest at L.
a = @(v, u) (-(P.a1 + P.m*P.a2)*v - P.a3*v.*abs(v) + P.b*u)/P.m;
vv = linspace(0, P.Vref, 200001);
ta = trapz(vv, 1./max(a(vv, 1), eps));   xa = trapz(vv, vv./max(a(vv, 1), eps));
tb = trapz(vv, 1./-a(vv, -1));            xb = trapz(vv, vv./-a(vv, -1));
tMin = ta + tb + (P.L - xa - xb)/P.Vref;
end
