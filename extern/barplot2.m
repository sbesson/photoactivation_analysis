%BARPLOT2 Bar plot grouping multiple sets/categories of data, with error bars.
%
% Inputs:   prm : cell array of matrices that contain the box properties:
%                 row 1: height
%                 row 2: optional, error bars
%     errorbars : cell array colors, each containing a Nx3 matrix. Colors cycle through matrix.
%
%       Options : see function content
%
% Examples: 
%
% 1) Simple bar plot
% figure; barplot2(rand(1,6), 0.1*rand(1,6), 'BarWidth', 0.8, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XTickLabel', arrayfun(@(k) ['S' num2str(k)], 1:6, 'UniformOutput', false),...
%     'Angle', 0, 'YLim', [0 1]);
%
% 2) Multiple groups
% figure; barplot2(rand(6,3), 0.1*rand(6,3), 'BarWidth', 1, 'GroupDistance', 1, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XTickLabel', arrayfun(@(k) ['group ' num2str(k)], 1:6, 'UniformOutput', false),...
%     'Angle', 45, 'YLim', [0 1.2]);
%
% 3) Multiple groups, asymmetric error bars
% figure;
% eb = 0.5*ones(2,4);
% barplot2([-2 1 -3 2; -2 1 -4 2], eb, eb/2, 'ErrorBarPosition', 'both',...
%     'BarWidth', 1, 'GroupDistance', 1, 'XLabel', 'x label', 'YLabel', 'y label');
%
%
% Note: this function uses patch() since colors can't be controlled with bar()
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

% Francois Aguet, 18 March 2011 (Last modified: 1 Mar 2013)

function he = barplot2(prm, varargin)

ng = size(prm,1); % #groups
nb = size(prm,2); % #bars in each group

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('prm');
ip.addOptional('errorbars', [], @(x) isempty(x) || all(size(x)==size(prm)));
ip.addOptional('BottomErrorbars', [], @(x) isempty(x) || all(size(x)==size(prm)));
ip.addOptional('Annotations', [], @(x) isempty(x) || size(x,2)==2);
ip.addParamValue('FaceColor', jet(nb), @(x) size(x,1)==1 || size(x,1)==nb || size(x,1)==ng);
ip.addParamValue('EdgeColor', []);
ip.addParamValue('BorderWidth', [], @isscalar); 
ip.addParamValue('XLabel', ' ', @ischar);
ip.addParamValue('YLabel', ' ', @ischar);
ip.addParamValue('YLim', [], @(x) numel(x)==2);
ip.addParamValue('YTick', []);
ip.addParamValue('XTickLabel', [], @(x) isempty(x) || (iscell(x) && (numel(x)==sum(nb) || numel(x)==ng)));
ip.addParamValue('BarWidth', 0.8, @isscalar);
ip.addParamValue('GroupDistance', 0.8, @isscalar);
ip.addParamValue('LineWidth', 1, @isscalar);
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('ErrorBarPosition', 'top',  @(x) strcmpi(x, 'top') | strcmpi(x, 'both'));
ip.addParamValue('ErrorBarWidth', 0.2, @(x) 0<x && x<=1);
ip.addParamValue('Handle', gca, @ishandle);
ip.addParamValue('Interpreter', 'tex', @(x) any(strcmpi(x, {'tex', 'latex', 'none'})));
ip.addParamValue('X', [], @(x) numel(x)==ng); % cell array of x-coordinates (groups only)
ip.addParamValue('AdjustFigure', true, @islogical);
ip.addParamValue('ErrorbarColor', []); % not yet implemented
ip.parse(prm, varargin{:});

topErrorbars = ip.Results.errorbars;
if isempty(ip.Results.BottomErrorbars)
    bottomErrorbars = topErrorbars;
else
    bottomErrorbars = ip.Results.BottomErrorbars;
end
faceColor = ip.Results.FaceColor;
if size(faceColor,1)==1
    faceColor = repmat(faceColor, [nb 1]);
end
nc = size(faceColor,1);

edgeColor = ip.Results.EdgeColor;
if size(edgeColor,1)==1
    edgeColor = repmat(edgeColor, [nb 1]);
elseif isempty(edgeColor)
    edgeColor = zeros(size(faceColor));
end

ha = ip.Results.Handle;
bw = ip.Results.BarWidth;
dg = ip.Results.GroupDistance; % distance between groups, in bar widths

% x-coords for groups
xa = cell(1,ng);
if isempty(ip.Results.X)
    xa{1} = 1:nb;
    for g = 2:ng
        xa{g} = xa{1} + xa{g-1}(end) + dg;
    end
else
    dx = min(diff(ip.Results.X));
    for g = 1:ng
        w = (nb-1)/2;
        xa{g} = ip.Results.X(g) + (g-1)*dg + (-w:w)*dx/nb;
    end
end
    
if isempty(ip.Results.BorderWidth)
    if ng>1
        border = bw/2+dg/2;
    else
        border = 1-bw/2;
    end
else
    border = ip.Results.BorderWidth;
end


hold on;
topval = zeros(1,ng*nb);
for g = 1:ng

    height = prm(g,:);
    
    % errorbars, if top only
    if ~isempty(topErrorbars)% && strcmpi(ip.Results.ErrorBarPosition, 'top')
        posIdx = height>=0;
        if sum(posIdx)>0
            he = errorbar(xa{g}(posIdx), height(posIdx), zeros(1,sum(posIdx)), topErrorbars(g,posIdx),...
                'k', 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
            setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'top');
            topval((g-1)*nb+find(posIdx)) = height(posIdx)+topErrorbars(g,posIdx);
        end
        posIdx = height<0; % negative values
        if sum(posIdx)>0
            he = errorbar(xa{g}(posIdx), height(posIdx), bottomErrorbars(g,posIdx), zeros(1,sum(posIdx)),...
                'k', 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
            setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'bottom');
        end
        %topval((g-1)*nb+find(posIdx)
    end
        
    % bars
    lb = xa{g} - bw/2;
    rb = xa{g} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [height; height; zeros(1,nb); zeros(1,nb); height; height];
    
    for b = 1:nb
        if nc==nb
            ci = b;
        else
            ci = g;
        end
        patch(xv(:,b), yv(:,b), faceColor(ci,:), 'EdgeColor', edgeColor(ci,:),...
            'LineWidth', ip.Results.LineWidth);
    end
   
    % errorbars, if two-sided
    if ~isempty(bottomErrorbars) && strcmpi(ip.Results.ErrorBarPosition, 'both')
        posIdx = height>=0;
        if sum(posIdx)>0
            he = errorbar(xa{g}(posIdx), height(posIdx), bottomErrorbars(g,posIdx), zeros(1,sum(posIdx)),...
                'k', 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
            setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'bottom');
        end
        posIdx = height<0; % negative values
        if sum(posIdx)>0
            he = errorbar(xa{g}(posIdx), height(posIdx), zeros(1,sum(posIdx)), topErrorbars(g,posIdx),...
                'k', 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
            setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'top');
        end
    end
end

% if there are negative values, plot axis
if min(prm(:)) < 0
    xAll = [xa{:}];
    plot([xAll(1)-border xAll(end)+border], [0 0], 'k-');
end

if ng>1
    XTick = arrayfun(@(k) (xa{k}(1) + xa{k}(end))/2, 1:ng);
else
    XTick = [xa{:}];
end

% position of the bars
xa = [xa{:}];

XTickLabel = ip.Results.XTickLabel;
if isempty(XTickLabel)
    if ng>1
        XTickLabel = 1:ng;
    else
        XTickLabel = 1:nb;
    end
end

XLim = [xa(1)-border xa(end)+border];
set(ha, 'XTick', XTick, 'XTickLabel', XTickLabel, 'XLim', XLim);
if ~isempty(ip.Results.YLim)
    YLim = ip.Results.YLim;
    set(ha, 'YLim', YLim);
else
    YLim = get(gca, 'YLim');
end
if ~isempty(ip.Results.YTick)
    set(ha, 'YTick', ip.Results.YTick);
end

% add annotation links (for significance etc.)
av = ip.Results.Annotations;
if ~isempty(av) && ~isempty(topErrorbars)
    pos = get(gca, 'Position');
    dy = 0.25/diff(XLim)*diff(YLim)/pos(4)*pos(3);
    maxposCount = zeros(numel(prm),1);
    for k = 1:size(av,1)
        y0 = max(topval(av(k,1):av(k,2)));
        maxpos = find(topval==y0, 1, 'first');
        maxposCount(maxpos) = maxposCount(maxpos)+1;
        plot(xa(av(k,[1 1 2 2])), y0+dy+1.75*dy*(maxposCount(maxpos)-1)+[0 dy dy 0], 'k',...
            'LineWidth', 0.75*ip.Results.LineWidth);
    end
end

% x labels
if ip.Results.Angle~=0 && ~isempty(ip.Results.XTickLabel)
    rotateXTickLabels(ha, 'Angle', ip.Results.Angle, 'Interpreter', ip.Results.Interpreter,...
        'AdjustFigure', ip.Results.AdjustFigure);
end
