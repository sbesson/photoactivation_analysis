function plotIndividualSeries(times, data, varargin)

% define small and large fonts for graphical output
default_tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
default_sfont = {'FontName', 'Helvetica', 'FontSize', 18};
default_lfont = {'FontName', 'Helvetica', 'FontSize', 22};

% Input checl
ip = inputParser;
ip.addRequired('times',@isvector);
ip.addRequired('data',@isnumeric);
ip.addOptional('events',[], @isvector);
ip.addParamValue('tfont', default_tfont, @iscell);
ip.addParamValue('sfont', default_sfont, @iscell);
ip.addParamValue('lfont', default_lfont, @iscell);
ip.parse(times, data, varargin{:});

% Plot indidividual events
plot(times, data);
nanData = ~all(isnan(data),2);
imin = max(1, find(nanData, 1, 'first'));
imax = min(find(nanData, 1, 'last'), size(data, 1));
xlim([times(imin) times(imax)]);
set(gca, 'LineWidth', 1.5, ip.Results.sfont{:}, 'Layer', 'top');
xlabel('Time (min)', ip.Results.lfont{:});
ylabel('Spindle length (\mum)', ip.Results.lfont{:})

if isempty(ip.Results.events), return; end
events = ip.Results.events;
nEvents = numel(events);
markers = 'osd';
hold on;
for i = 1 : nEvents
    plot(events(i).times, events(i).values, 'k',...
        'Marker', markers(i), 'MarkerFaceColor', 'k');
end
legend({events.name},ip.Results.tfont{:})
;
