function [f, P] = lyapunov_margin(A1, A2)
%LYAPUNOV_MARGIN Common quadratic Lyapunov function for two discrete-time matrices.
%   [f, P] = LYAPUNOV_MARGIN(A1, A2) returns
%       f = min over P >= 0, trace(P) = 1, of  max_i lambda_max(Ai'*P*Ai - P)
%   and the minimizing P. A common P > 0 with Ai'*P*Ai - P < 0 (i = 1, 2)
%   exists exactly when f < 0 (the conditions are homogeneous in P, so the
%   trace can be fixed). The function inside the min is convex in P, so a
%   grid search followed by a local refinement finds the global minimum.
%   Plain MATLAB, no LMI solver needed.
%
%   P = [p, q; q, 1 - p] is parameterized as p = (1 + tanh(s1))/2,
%   q = tanh(s2)*sqrt(p*(1 - p)), which covers every P >= 0 with trace 1.

obj = @(s) worst_eig(s, A1, A2);
best = inf; s0 = [0 0];
for s1 = linspace(-4, 4, 81)
    for s2 = linspace(-4, 4, 81)
        val = obj([s1 s2]);
        if val < best, best = val; s0 = [s1 s2]; end
    end
end
opts = optimset('TolX', 1e-12, 'TolFun', 1e-14, 'MaxFunEvals', 1e5, 'MaxIter', 1e5);
s = fminsearch(obj, s0, opts);
f = obj(s);
P = make_P(s);
end

function val = worst_eig(s, A1, A2)
P = make_P(s);
val = max(max(eig(A1'*P*A1 - P)), max(eig(A2'*P*A2 - P)));
end

function P = make_P(s)
p = (1 + tanh(s(1)))/2;
q = tanh(s(2))*sqrt(p*(1 - p));
P = [p q; q 1 - p];
end
