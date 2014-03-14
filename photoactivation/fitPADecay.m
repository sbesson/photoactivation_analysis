function [bFit, r2, se] = fitPADecay(times, dI, varargin)

% Parse input range
ip = inputParser;
ip.addOptional('maxTime', max(times), @isscalar);
ip.parse(varargin{:});
range = times <= ip.Results.maxTime;

%Double-exponential function for fitting
turnoverFn = @(b,x) b(1) .* exp(-x/b(2)) + b(3) .* exp(-x/b(4)) + (1-b(1)-b(3));
diffFn = @(p) turnoverFn(p, times(range)) - dI(range);


% Initial parameters
bInit = [.8 20 .2 100];
lb = [0 1 0 40];
ub = [1.2 40 1.2 600];

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
% ;
% [bFit,r,J] = nlinfit(times(range), dI(range), turnoverFn,...
%     bInit);
% resnorm=sum(r.^2);
r2 = 1 - resnorm / norm(dI(range)-mean(dI(range)))^2;


alpha = 0.05;
ci = nlparci(bFit,r,'Jacobian',J,'alpha',alpha);
t = tinv(1-alpha/2,length(dI(range))-length(bFit));
se = (ci(:,2)-ci(:,1)) ./ (2*t);  % Standard Error
