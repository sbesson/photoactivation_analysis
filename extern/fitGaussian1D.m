function [prmVect, prmStd, C, res, J] = fitGaussian1D(x, data, prmVect, mode, mask)
%[prmVect, G] = fitGaussian2D(data, p, mode)
%
% Input: data: 2-D image array
%        p      : [xp yp A sigma c] initial and fixed parameter values
%        {mode} : specifies which parameters to estimate; any combination of 'xAsc'
%        {mask} : elements set to 1 are not included in optimization
%        {xa}   : x-axis
%
% Data is assumed to contain a single Gaussian
%
% Francois Aguet, March 30 2011

%
% Copyright (C) 2014 LCCB 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
%

if nargin<3
    mode = 'xAsc';
end
if nargin<5
    mask = [];
end

opts = optimset('Jacobian', 'on', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);


estIdx = false(1,4); % [x A s c]
estIdx(regexpi('xAsc', ['[' mode ']'])) = true;
lb = [x(1) 0 0 0];
ub = [x(end) Inf Inf Inf];


[p, resnorm, res, ~, ~, ~, J] = lsqnonlin(@costGaussian, prmVect(estIdx), lb(estIdx), ub(estIdx), opts, data, x, prmVect, estIdx, mask);
%[p, resnorm, res, ~, ~, ~, J] = lsqnonlin(@costGaussian, prmVect(estIdx), [], [], opts, data, x, y, prmVect, estIdx, mask);
prmVect(estIdx) = p;


sigma2 = resnorm / (numel(data) - sum(estIdx) - 1);
J = full(J);
C = inv(J'*J);
prmStd = sqrt(sigma2*diag(C))';


function [v, J] = costGaussian(p, data, x, prmVect, estIdx, mask)
prmVect(estIdx) = p;

[g J] = gaussian1D(x, prmVect);
J(:,estIdx==false) = []; % remove unneeded Jacobian components
v = g - data;
maskIdx = mask==1;
v(maskIdx) = [];
J(maskIdx, :) = [];


function [g J] = gaussian1D(x, prmVect)

tmp = num2cell(prmVect);
[xp A s c] = deal(tmp{:});

r2 = (x-xp).^2;

g_dA = exp(-r2/(2*s^2));
g = A*g_dA;

g_dxp = (x-xp)./s^2.*g;
g_ds = r2/s^3.*g;
g = g + c;

N = numel(x);
g_dc = ones(N,1);
J = [reshape(g_dxp, [N 1]) reshape(g_dA, [N 1]) reshape(g_ds, [N 1]) g_dc];
