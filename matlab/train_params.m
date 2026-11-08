function P = train_params(varargin)
%TRAIN_PARAMS Parameters of the electric train and of the trip (Theme 1).
%   P = TRAIN_PARAMS() returns the nominal values.
%   P = TRAIN_PARAMS('m', 110, 'grade', 0.02, ...) overrides fields by name.
%
%   Model:  m*dv/dt = -(a1 + m*a2)*v - a3*v*|v| + b*u - m*g*sin(theta),
%           dx/dt = v,  |u| <= 1,  tan(theta) = grade (+ uphill, - downhill).
%   Units: m, s, kg, N.

P.m  = 100;     % mass (nominal)
P.a1 = 1;       % friction parameters
P.a2 = 0.05;
P.a3 = 0.1;
P.b  = 200;     % motor gain
P.g  = 9.81;
P.grade = 0;    % slope of the track, rise/run (0.02 = 2 %)

P.L    = 30000; % distance between the two stations (m)
P.Vref = 20;    % reference speed between the stations (m/s, 72 km/h)
P.Ts   = 0.1;   % controller sample time (s); the plant is integrated with RK4 at the same step
P.Tend = 2400;  % simulated time (s)

if mod(numel(varargin), 2) ~= 0
    error('train_params:args', 'Parameters must be given as name-value pairs.');
end
for k = 1:2:numel(varargin)
    if ~isfield(P, varargin{k})
        error('train_params:unknown', 'Unknown parameter "%s".', varargin{k});
    end
    P.(varargin{k}) = varargin{k+1};
end
end
