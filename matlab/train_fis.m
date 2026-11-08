function [fis, C] = train_fis(varargin)
%TRAIN_FIS Takagi-Sugeno fuzzy controller for the train (Theme 1).
%   [fis, C] = TRAIN_FIS() builds the controller with the default settings.
%   [fis, C] = TRAIN_FIS('adec', 0.5, 'nNodes', 2, ...) changes them.
%
%   Inputs:  d = distance to the next station (m),  v = speed (m/s).
%   The braking zone is covered by fuzzy sets D_1 ... D_n centred at the
%   distances d_i where a train braking at the constant rate adec still has
%   the speed v_i:  v_i = i*Vref/n,  d_i = v_i^2/(2*adec).  Rules:
%     Past  IF d is Past THEN u = -uff(v_1) + K*(-v_1 - v) (d < 0: mirror of D_1, back to the station)
%     Stop  IF d is Stop THEN u = -K*v                     (at the station, v_0 = 0)
%     D_i   IF d is D_i  THEN u = uff(v_i) + K*(v_i - v)   (i = 1 ... n-1)
%     Far   IF d is Far  THEN u = uff(Vref) + Kv*(Vref - v)
%   uff(v) = ((a1 + m*a2)*v + a3*v^2)/b is the input that holds the speed v
%   on level track with the nominal mass. With triangular sets that add up
%   to 1 and the same gain in every rule (Kv = K), the weighted average is
%   u = uff + K*(vhat(d) - v), where the target speed vhat(d) is the
%   piecewise-linear interpolation of the centres (d_i, v_i), close to
%   min(Vref, sqrt(2*adec*d)). The caller clips u to [-1, 1].
%
%   Settings (name-value):
%     Vref        reference speed between the stations (m/s)   default 20
%     adec        braking rate of the speed profile (m/s^2)     default 1
%     nNodes      number of speed steps from 0 to Vref          default 5
%     K           speed gain of the braking-zone rules (s/m)    default 1.25
%     Kv          speed gain of the cruise rule (s/m)           default 1.25 (= K)
%     ff          use the feed-forward terms uff                default true
%     shape       'tri' (sets add up to 1) or 'gauss'           default 'tri'
%     defuzz      'wtaver' or 'wtsum'                           default 'wtaver'

C.Vref = 20; C.adec = 1; C.nNodes = 5; C.K = 1.25; C.Kv = 1.25; C.ff = true;
C.shape = 'tri'; C.defuzz = 'wtaver';
for k = 1:2:numel(varargin)
    if ~isfield(C, varargin{k})
        error('train_fis:unknown', 'Unknown setting "%s".', varargin{k});
    end
    C.(varargin{k}) = varargin{k+1};
end

P = train_params();
uff = @(v) C.ff*((P.a1 + P.m*P.a2)*v + P.a3*v.^2)/P.b;
n = C.nNodes;
C.vNodes = (0:n)*C.Vref/n;                 % node speeds, the last one is Vref
C.dNodes = C.vNodes.^2/(2*C.adec);         % node distances, the last one starts Far
d = C.dNodes; Dmax = 1e6;

fis = sugfis('Name', 'train', 'DefuzzificationMethod', C.defuzz);
fis = addInput(fis, [-Dmax Dmax], 'Name', 'd');
fis = addInput(fis, [-100 100], 'Name', 'v');
fis = addMF(fis, 'v', 'trapmf', [-200 -100 100 200], 'Name', 'Any');
fis = addOutput(fis, [-1 1], 'Name', 'u');
names = [{'Past', 'Stop'}, arrayfun(@(i) sprintf('D%d', i), 1:n-1, 'UniformOutput', false), {'Far'}];
switch C.shape
    case 'tri'
        fis = addMF(fis, 'd', 'trapmf', [-2*Dmax -Dmax -d(2) 0], 'Name', 'Past');
        fis = addMF(fis, 'd', 'trimf', [-d(2) 0 d(2)], 'Name', 'Stop');
        for i = 1:n-1
            fis = addMF(fis, 'd', 'trimf', [d(i) d(i+1) d(i+2)], 'Name', names{i+2});
        end
        fis = addMF(fis, 'd', 'trapmf', [d(n) d(n+1) Dmax 2*Dmax], 'Name', 'Far');
    case 'gauss'
        % Gaussian sets at the same centres, with a width of half the gap to
        % the nearer neighbour. They never reach 0 and do not add up to 1.
        gap = diff(d);
        s = [gap(1), min(gap(1:end-1), gap(2:end)), gap(end)]/2;   % widths at d_0 ... d_n
        fis = addMF(fis, 'd', 'gauss2mf', [1 -2*Dmax s(1) -d(2)], 'Name', 'Past');
        fis = addMF(fis, 'd', 'gaussmf', [s(1) 0], 'Name', 'Stop');
        for i = 1:n-1
            fis = addMF(fis, 'd', 'gaussmf', [s(i+1) d(i+1)], 'Name', names{i+2});
        end
        fis = addMF(fis, 'd', 'gauss2mf', [s(n+1) d(n+1) 1 2*Dmax], 'Name', 'Far');
    otherwise
        error('train_fis:shape', 'Unknown shape "%s".', C.shape);
end
% Consequents u = c_d*d + c_v*v + c_0.
fis = addMF(fis, 'u', 'linear', [0, -C.K, -(uff(C.vNodes(2)) + C.K*C.vNodes(2))], 'Name', 'Past');
fis = addMF(fis, 'u', 'linear', [0, -C.K, 0], 'Name', 'Stop');
for i = 1:n-1
    fis = addMF(fis, 'u', 'linear', [0, -C.K, uff(C.vNodes(i+1)) + C.K*C.vNodes(i+1)], 'Name', names{i+2});
end
fis = addMF(fis, 'u', 'linear', [0, -C.Kv, uff(C.Vref) + C.Kv*C.Vref], 'Name', 'Far');
rules = strings(1, numel(names));
for i = 1:numel(names)
    rules(i) = sprintf('d==%s => u=%s', names{i}, names{i});
end
fis = addRule(fis, rules);
C.names = names;
C.uffRef = uff(C.Vref);
% Linearization at the station (0 < d < d_1, small v): u = Kp*d - K*v with
% Kp = (uff(v_1) + K*v_1)/d_1, so m*d'' + (a1 + m*a2 + b*K)*d' + b*Kp*d = 0.
C.Kp = (uff(C.vNodes(2)) + C.K*C.vNodes(2))/d(2);
C.zeta = @(m) (P.a1 + m*P.a2 + P.b*C.K)/(2*sqrt(m*P.b*C.Kp));
end
