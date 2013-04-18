function [events, log] = detectEvents(data,times, minSpindleLength)
% Detects multiple events in time-series
%
% SYNOPSIS [events, log] = detectEvents(data, times, minSpindleLength)
%
% This function detects individual events in timeseries.
% 1 - the onset, i.e. where the value exceeds a given threshold.
% 2 - the maximum
%
% INPUT
%       data - a nTimepoints x nSamples matrix containing the time-series
%
%       times - a nTimepoints x 1 matrix containing the time points
%
%       minSpindleLengh - the minimum spindle length to detect events and
%       align the time series
%
% OUTPUT
%      events - an array of structures containing the events. Each element
%      of the array contains the following field
%          
%          name: Name of the detection event
%          detectionFunc: Function applied to detect the event
%          index: index of the event for each sample
%          values: value of the events for each sample
%          times: time where the event occurs for each sample
%
%      log - a string containing the log of the event detection function

% Sebastien Besson Dec 2012 (last modified Apr 2013)

% Input check
ip = inputParser;
ip.addRequired('data',@isnumeric);
nTimePoints = size(data, 1);
nSamples = size(data, 2);
ip.addRequired('times',@(x) isvector(x) && numel(x) == nTimePoints);
ip.addRequired('minSpindleLength', @isscalar);
ip.parse(data, times, minSpindleLength);

%% Initialize output structure array

events = deal(struct('name', '', 'detectionFunc', [],...
    'index', NaN(1, nSamples), 'values', NaN(1, nSamples), 'times', NaN(1, nSamples)));

%% Elongation onset detection
events(1).name = 'Spindle elongation onset';
events(1).detectionFunc = @(x) max(1, find(x > minSpindleLength, 1, 'first')-1);
log = sprintf('Finding last events where value exceeds %g\n', minSpindleLength);

% Detect events for each time-series (values above threshold)
for i = 1:nSamples
    index = events(1).detectionFunc(data(:,i));
    if ~isempty(index)
        events(1).index(i) = index;
        events(1).values(i) = data(events(1).index(i), i);
    end
end

%% Maximum elongation
events(2).name = 'Maximum spindle elongation';
events(2).detectionFunc = @(x) max(x);
log = [log sprintf('Finding maximum events\n')];

% Detect events for each time-series (values above threshold)
for i = 1:nSamples
    [events(2).values(i), events(2).index(i)] = events(2).detectionFunc(data(:,i));
end

%% Maximum elongation
events(3).name = 'Spindle elongation at 10 min';
log = [log sprintf('Finding maximum events\n')];

% Detect events for each time-series (values above threshold)
events(3).index = min(events(1).index + 2, nTimePoints);
events(3).values = NaN(1, nSamples);
for i = find(~isnan(events(3).index))
    events(3).values(i) =  data(events(3).index(i), i);
end

%% Event times
for i = 1:numel(events)
    validEvents = ~isnan(events(i).index);
    events(i).times(validEvents) = times(events(i).index(validEvents));
    events(i).times(~validEvents) = NaN;
    % Read event times
    log = [log sprintf([events(i).name ' mean time: %g +/- %g\n'], ...
        mean(events(i).times), std(events(i).times,[],2))];
end