function [bFit, r2, ci] = fitPADecay(times, dI, varargin)

% Parse input range
ip = inputParser;
ip.addOptional('maxTime', max(times), @isscalar);
ip.parse(varargin{:});
range = times <= ip.Results.maxTime;

%Double-exponential function for fitting
turnoverFn = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x)) + (1-b(1)-b(3));
diffFn = @(p) turnoverFn(p, times(range)) - dI(range);


% Initial parameters
bInit = [.8 -.1 .2 -0.01];
lb = [0 -Inf 0 -Inf];
ub = [1.2 0 1.2 0];

% Define fitting options
opts = optimset('MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

% Use non-linear least squares
[bFit, resnorm,r, ~,~,~, J] = lsqnonlin(diffFn, bInit, lb, ub, opts);
% [bFit, resnorm, r, ~, ~, J] =  lsqcurvefit(...
%     turnoverFn, bInit, times(range), dI(range));
% resnorm=sum(r.^2);
r2 = 1 - resnorm / norm(dI-mean(dI))^2;


ci = nlparci(bFit,r,'Jacobian',J);

