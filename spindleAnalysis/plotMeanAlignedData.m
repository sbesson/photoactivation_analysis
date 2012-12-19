function plotMeanAlignedData(s, varargin)

% define small and large fonts for graphical output
default_tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
default_sfont = {'FontName', 'Helvetica', 'FontSize', 18};
default_lfont = {'FontName', 'Helvetica', 'FontSize', 22};

% Input checl
ip = inputParser;
ip.addRequired('s',@isstruct);
ip.addOptional('nPointsMin', 5, @isscalar);
ip.addParamValue('tfont', default_tfont, @iscell);
ip.addParamValue('sfont', default_sfont, @iscell);
ip.addParamValue('lfont', default_lfont, @iscell);
ip.parse(s, varargin{:});

% Read number of elements to display
nElements = numel(s);
colors = hsv(nElements);
tmin = Inf;
tmax = -Inf;
time_max = max([s.alignedTimes]);

for i = 1 : nElements
    
    % Calculate mean and std values of aligned data
    alignedmean = nanmean(s(i).alignedData,2);
    alignedstd = nanstd(s(i).alignedData, [], 2);
    nPoints  =sum(~isnan(s(i).alignedData),2);
    alignedste = alignedstd./sqrt(nPoints);
    
    % Filter series with too few points
    exclude = nPoints < ip.Results.nPointsMin;
    alignedmean(exclude) = NaN;
    alignedste(exclude) = NaN;
    
    % Calculate mininimum and maximum times
    imin = max(1,find(~isnan(alignedmean),1,'first'));
    tmin = min(tmin, s(i).alignedTimes(imin));
    imax = min(find(~isnan(alignedmean(imin:end)),1,'last'), time_max);
    tmax = max(tmax, s(i).alignedTimes(imin+imax));

    hold on;
    errorbar(s(i).alignedTimes', alignedmean, alignedste,...
        '-','Color', colors(i,:), 'LineWidth', 1);
end

% Set axes properties
set(gca, 'LineWidth', 1.5, ip.Results.sfont{:}, 'Layer', 'top');
xlabel('Time (min)', ip.Results.lfont{:});
ylabel('Spindle length (\mum)', ip.Results.lfont{:});
xlim([tmin tmax]);
legend({s.name},'Location', 'NorthWest','Interpreter','None')
set(gca,'LooseInset',get(gca,'TightInset'));
box on