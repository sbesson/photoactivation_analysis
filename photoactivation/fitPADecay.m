function [bFit, r2] = fitPADecay(times, dI, varargin)

% Parse input range
ip = inputParser;
ip.addOptional('maxTime', max(times), @isscalar);
ip.parse(varargin{:});
range = times <= ip.Results.maxTime;

%Double-exponential function for fitting
turnoverFn = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x)) +b(5);
diffFn = @(p) turnoverFn(p, times(range)) - dI(range);


% Initial parameters
bInit = [.8 -.1 .2 -0.01 0];
lb = [0 -Inf 0 -Inf 0];
ub = [1.2 0 1.2 0 .2];

% Define fitting options
opts = optimset('MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

% Use non-linear least squares
[bFit, resnorm] = lsqnonlin(diffFn, bInit, lb, ub, opts);
r2 = 1 - resnorm / norm(dI-mean(dI))^2;

