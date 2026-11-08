function R = simulate_train(fis, P)
%SIMULATE_TRAIN Closed-loop simulation of the train with a fuzzy controller.
%   R = SIMULATE_TRAIN(fis, P) starts at rest at x = 0 and drives to x = P.L.
%   The controller is sampled every P.Ts seconds (zero-order hold) and its
%   output is clipped to [-1, 1]. The plant is integrated with one RK4 step
%   per sample. R holds t, x, v, u and the summary measures:
%     tArrive   first time after which |L - x| <= 0.5 m and |v| <= 0.05 m/s
%     overshoot max(x) - L (m); a positive value means the train passed the mark
%     stopErr   final position minus L (m)
%     vCruise   mean speed while 5 km < x < L - 5 km
%     vMax      maximum speed;  aMax, aMin  extreme accelerations (m/s^2)
%     uMax      max |u|;  duMax  largest change of u between two samples

opt = evalfisOptions('OutOfRangeInputValueMessage', 'none', ...
    'NoRuleFiredMessage', 'none', 'EmptyOutputFuzzySetMessage', 'none');
N = round(P.Tend/P.Ts);
t = (0:N)*P.Ts;
x = zeros(1, N+1); v = x; u = x; acc = x;
f = @(s, uk) [s(2); plant_acc(s(2), uk, P)];
s = [0; 0];
for k = 1:N+1
    uk = evalfis(fis, [P.L - s(1), s(2)], opt);
    uk = min(max(uk, -1), 1);
    x(k) = s(1); v(k) = s(2); u(k) = uk; acc(k) = plant_acc(s(2), uk, P);
    if k <= N
        h = P.Ts;
        k1 = f(s, uk); k2 = f(s + h/2*k1, uk); k3 = f(s + h/2*k2, uk); k4 = f(s + h*k3, uk);
        s = s + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end

R.t = t; R.x = x; R.v = v; R.u = u; R.acc = acc;
ok = abs(P.L - x) <= 0.5 & abs(v) <= 0.05;
last = find(~ok, 1, 'last');
if isempty(last), R.tArrive = 0;
elseif last == numel(t), R.tArrive = NaN;          % never settles
else, R.tArrive = t(last + 1);
end
R.overshoot = max(x) - P.L;
R.stopErr = x(end) - P.L;
mid = x > 5000 & x < P.L - 5000;
R.vCruise = mean(v(mid));
R.vMax = max(v);
R.aMax = max(acc); R.aMin = min(acc);
R.uMax = max(abs(u));
R.duMax = max(abs(diff(u)));                % largest change of u in one sample
end

function a = plant_acc(v, u, P)
a = (-(P.a1 + P.m*P.a2)*v - P.a3*v*abs(v) + P.b*u)/P.m - P.g*sin(atan(P.grade));
end
