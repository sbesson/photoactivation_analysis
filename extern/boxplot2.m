%BOXPLOT2 Box plot grouping multiple sets/categories of data, with error bars and SEM.
%
% Inputs:   prm : matrix or cell array of matrices that contain the box properties:
%                 row 1: mean or median
%                 row 2: optional, SEM
%                 row 2/3: 25th percentile, bottom of box
%                 row 3/4: 75th percentile, top of box
%                 row 4/5: optional, bottom whisker
%                 row 5/6: optional, top whisker
% Options:
%  
%     FaceColor : Nx3 matrix of colors, where N is the number of bars or groups
%     EdgeColor : "
%       xLabels : cell array of strings, labels for each bar
%        yLabel : string, y-axis label
%
% Examples:
%
% 1) Simple box plot
% prm = [3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5];
% figure; boxplot2(prm, 'BarWidth', 0.8, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XTickLabel', arrayfun(@(k) ['S' num2str(k)], 1:2, 'UniformOutput', false),...
%     'Angle', 0);
% 
% 
% 2) Multiple groups
% prm = {[3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5],...
%     [3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5]};
% figure; boxplot2(prm, 'BarWidth', 0.8, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XTickLabel', arrayfun(@(k) ['Group label ' num2str(k)], 1:2, 'UniformOutput', false),...
%     'Angle', 45, 'FaceColor', [1 0.5 0.5; 0.5 1 0.5], 'EdgeColor', [0.8 0 0; 0 0.8 0]);

% Francois Aguet, 22 Feb 2011 (Last modified: 08/15/2012)

function h = boxplot2(prm, varargin)

if isnumeric(prm)
    prm = {prm};
end

nbin = size(prm{1},2); % # bins
nd = numel(prm); % # data sets

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('prm');
ip.addOptional('Annotations', [], @(x) isempty(x) || size(x,2)==2);
ip.addParamValue('FaceColor', jet(max(nd)), @(x) size(x,1)==1 || size(x,1)==nd || size(x,1)==nbin);
ip.addParamValue('EdgeColor', zeros(1,3));
ip.addParamValue('GroupDistance', 0.5, @isscalar);
ip.addParamValue('BorderWidth', [], @isscalar); 
ip.addParamValue('XLabel', [], @ischar);
ip.addParamValue('XTickLabel', arrayfun(@(k) num2str(k), 1:sum(nbin), 'UniformOutput', false), @(x) iscell(x) && (numel(x)==sum(nd)||numel(x)==nbin));
ip.addParamValue('YLabel', ' ', @ischar);
ip.addParamValue('YLim', [], @(x) numel(x)==2);
ip.addParamValue('BarWidth', 0.8, @isscalar);
ip.addParamValue('LineWidth', 1, @isscalar);
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('ErrorBarWidth', 0.2, @(x) 0<x && x<=1);
ip.addParamValue('Handle', gca, @ishandle);
ip.addParamValue('Interpreter', 'tex', @(x) any(strcmpi(x, {'tex', 'latex', 'none'})));
ip.addParamValue('X', [], @(x) numel(x)==nbin); % cell array of x-coordinates (groups only)
ip.addParamValue('AdjustFigure', true, @islogical);
ip.addParamValue('ErrorbarColor', []);
ip.addParamValue('PlotDensity', false, @islogical);
ip.parse(prm, varargin{:});

faceColor = ip.Results.FaceColor;
if ~iscell(faceColor) % otherwise assume correct format: cell(1,nd)(nbin,:)
    if size(faceColor,1)==1
        faceColor = repmat({repmat(faceColor, [nbin 1])}, [1 nd]);
    elseif size(faceColor,1)==nbin
        faceColor = repmat({faceColor}, [1 nd]);
    elseif size(faceColor,1)==nd
        faceColor = arrayfun(@(i) repmat(faceColor(i,:), [nbin 1]), 1:nd, 'UniformOutput', false);
    else
        error('Incorrect color format.');
    end
end

edgeColor = ip.Results.EdgeColor;
if ~iscell(edgeColor)
    if size(edgeColor,1)==1
        edgeColor = repmat({repmat(edgeColor, [nbin 1])}, [1 nd]);
    elseif size(edgeColor,1)==nbin
        edgeColor = repmat({edgeColor}, [1 nd]);
    elseif size(edgeColor,1)==nd
        edgeColor = arrayfun(@(i) repmat(edgeColor(i,:), [nbin 1]), 1:nd, 'UniformOutput', false);
    else
        error('Incorrect color format.');
    end
end

errorbarColor = ip.Results.ErrorbarColor;
if isempty(errorbarColor)
    errorbarColor = repmat({zeros(nbin,3)}, [1 nd]);
end

ha = ip.Results.Handle;
bw = ip.Results.BarWidth;
dg = ip.Results.GroupDistance; % distance between groups, in bar widths

% x-coords for groups
xa = cell(1,nbin);
if isempty(ip.Results.X)
    xa{1} = 1:nd;
    for k = 2:nbin
        xa{k} = xa{1} + xa{k-1}(end) + dg;
    end
else
    dx = min(diff(ip.Results.X));
    for k = 1:nbin
        w = (nd-1)/2;
        xa{k} = ip.Results.X(k) + (k-1)*dg + (-w:w)*dx/nd;
    end
end

if isempty(ip.Results.BorderWidth)
    border = 0.5+dg/2;
else
    border = ip.Results.BorderWidth;
end

if ~iscell(prm{1})
    plotSEM = mod(size(prm{1},1),2)==0;
else
    plotSEM = false;
end

hold on;
% handles
h = zeros(1,nd);
topval = zeros(1,nbin*nd);
for k = 1:nbin
    
    % concatenate values for group 'k'
    if iscell(prm{1})
        M = cellfun(@(i) [prctile(i{k}, [50 25 75]) min(i{k}) max(i{k})], prm, 'UniformOutput', false);
        M = vertcat(M{:})';
    else
        M = cellfun(@(i) i(:,k), prm, 'UniformOutput', false);
        M = [M{:}];
    end
        
    if plotSEM
        p25 = M(3,:);
        p75 = M(4,:);
    else
        p25 = M(2,:);
        p75 = M(3,:);
    end
    
    % whiskers (plot first to mask bar at '0')
    if plotSEM && size(M,1)>4
        w1 = M(5,:);
        w2 = M(6,:);
        plotWhiskers = 1;
    elseif size(M,1)>3
        w1 = M(4,:);
        w2 = M(5,:);
        plotWhiskers = 1;
    else
        plotWhiskers = 0;
    end
    
    mu = M(1,:);
    
    % the box
    lb = xa{k} - bw/2;
    rb = xa{k} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [p75; p75; p25; p25; p75; p75];
    topval((k-1)*nd+(1:nd)) = M(end,:);
    
    for b = 1:nd

        if ip.Results.PlotDensity
            [f,xi] = ksdensity(prm{b}{k});
            f = f/max(f)/2.5;
            fill([xa{k}(b)+f xa{k}(b)-f(end:-1:1)], [xi xi(end:-1:1)], faceColor{b}(k,:),...
                'EdgeColor', edgeColor{b}(k,:), 'HandleVisibility', 'off');
        end
        
        if plotWhiskers
            he = errorbar(xa{k}(b), p25(b), w1(b)-p25(b), 0, 'Color', errorbarColor{b}(k,:), 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
            setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'bottom');
            he = errorbar(xa{k}(b), p75(b), 0, w2(b)-p75(b), 'Color', errorbarColor{b}(k,:), 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
            setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'top');
        end
        
        hp = patch(xv(:,b), yv(:,b), faceColor{b}(k,:), 'EdgeColor', edgeColor{b}(k,:),...
            'LineWidth', ip.Results.LineWidth);
        if k==1
            h(b) = hp;
        end
        
        % mean/median line
        line([lb(b); rb(b)], [mu(b); mu(b)], 'Color', errorbarColor{b}(k,:), 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
        
        % replot border
        hp = patch(xv(:,b), yv(:,b), faceColor{b}(k,:), 'EdgeColor', edgeColor{b}(k,:),...
            'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
        set(hp, 'FaceColor', 'none');
    end
    
    
    % SEM
    if plotSEM
        sigma = M(2,:);
        he = errorbar(xa{k}, mu, sigma, 'k', 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
        setErrorbarStyle(he, 0.15);
    end
end
box off;

if numel(ip.Results.XTickLabel)==nbin
    XTick = arrayfun(@(k) (xa{k}(1) + xa{k}(end))/2, 1:nbin);
else
    XTick = [xa{:}];
end

% position of the bars
xa = [xa{:}];

XLim = [xa(1)-border xa(end)+border];
set(ha, 'XTick', XTick, 'XTickLabel', ip.Results.XTickLabel, 'XLim', XLim);
if ~isempty(ip.Results.YLim);
    YLim = ip.Results.YLim;
    set(ha, 'YLim', YLim);
else
    YLim = get(gca, 'YLim');
end

% add annotation links (for significance etc.)
av = ip.Results.Annotations;
if ~isempty(av)
    pos = get(gca, 'Position');
    dy = 0.25/diff(XLim)*diff(YLim)/pos(4)*pos(3);
    maxposCount = zeros(nbin,1);
    for k = 1:size(av,1)
        y0 = max(topval(av(k,1):av(k,2)));
        maxpos = find(topval==y0, 1, 'first');
        maxposCount(maxpos) = maxposCount(maxpos)+1;
        plot(xa(av(k,[1 1 2 2])), y0+dy+1.75*dy*(maxposCount(maxpos)-1)+[0 dy dy 0], 'k', 'LineWidth', 0.75);
    end
end
hold off;

% x labels
if ip.Results.Angle ~= 0
    rotateXTickLabels(ha, 'Angle', ip.Results.Angle, 'Interpreter', ip.Results.Interpreter,...
        'AdjustFigure', ip.Results.AdjustFigure);
end
