%ROTATEXTICKLABELS rotates the labels on the x-axis
%
% INPUTS:      ha : axis handle of the plot to be modified
%         'Angle' : rotation in degrees. Default: 45

% Francois Aguet, 22 Feb 2011
% Last modified: 08/15/2012

function ht = rotateXTickLabels(ha, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('ha', @ishandle);
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('Interpreter', 'tex', @(x) any(strcmpi(x, {'tex', 'latex', 'none'})));
ip.addParamValue('AdjustFigure', true, @islogical);
ip.parse(ha, varargin{:});

xa = get(ha, 'XTick');
xla = get(ha, 'XTickLabel');
if ischar(xla) % label is numerical
    xla = arrayfun(@(i) num2str(str2double(xla(i,:))), 1:size(xla,1), 'UniformOutput', false);
end
axPos = get(ha, 'Position');
set(ha, 'XTickLabel', [], 'Position', axPos);

fontName = get(ha, 'FontName');
fontSize = get(ha, 'FontSize');

XLim = get(ha, 'XLim');
YLim = get(ha, 'Ylim');
width = diff(XLim);
height = diff(YLim);


% get height of default text bounding box
h = text(0, 0, ' ', 'FontName', fontName, 'FontSize', fontSize);
extent = get(h, 'extent');
shift = extent(4)/height*width/axPos(3)*axPos(4) * sin(ip.Results.Angle*pi/180)/4;
delete(h);


ht = arrayfun(@(k) text(xa(k)-shift, YLim(1)-0.02*height, xla{k},...
    'FontName', fontName, 'FontSize', fontSize,...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
    'Rotation', ip.Results.Angle, 'Interpreter', ip.Results.Interpreter, 'Parent', ha),...
    1:length(xa));

% largest extent (relative to axes units)
extents = arrayfun(@(k) get(k, 'extent'), ht, 'UniformOutput', false);
extents = vertcat(extents{:});


lmargin = -min(extents(:,1))/width * axPos(3); % normalized units in fig. frame

hx = get(ha, 'XLabel');
maxHeight = max(extents(:,4));
if ~strcmpi(get(hx, 'String'), ' ')
    xlheight = get(hx, 'Extent');
    xlheight = xlheight(4);
else
    xlheight = 0;
end
bmargin = (maxHeight+xlheight)/height * axPos(4); % data units -> normalized

if ip.Results.AdjustFigure
    hfig = get(ha, 'Parent');
    fpos = get(hfig, 'Position');
    % expand figure window
    
    if lmargin > axPos(1)
        fpos(3) = fpos(3) + lmargin-axPos(1);
        axPos(1) = lmargin;
    end
    
    switch get(ha, 'Units')
        case 'normalized'
            cf = 1-axPos(2) + bmargin;
            axPos(4) = axPos(4) / cf;
            fpos(4) = fpos(4) * cf;
            axPos(2) = bmargin / cf;
        otherwise % assumes units are the same for figure and axes
            % extend figure, add bottom offest
            fpos(4) = fpos(4) - axPos(2) + bmargin;
            axPos(2) = bmargin;
    end
    
    set(hfig, 'Position', fpos, 'PaperPositionMode', 'auto');
    set(ha, 'Position', axPos);
end

% shift x-label
xPos = get(hx, 'Position');
xPos(2) = -maxHeight;
set(hx, 'Position', xPos, 'VerticalAlignment', 'middle');
