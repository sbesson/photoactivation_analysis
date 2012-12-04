function [x] = norminv(p, mu, sigma)
% This function mimicks the simple functionality of norminv


if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end

x0 = -sqrt(2).*erfcinv(2*p);
x = sigma.*x0 + mu;