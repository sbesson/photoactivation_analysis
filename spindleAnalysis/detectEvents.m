function [events, log] = detectEvents(data, minSpindleLength)
% Detects multiple events in time-series
%
% SYNOPSIS [events, log] = detectEvents(data, minSpindleLength)
%
% This function Then it detects the
% timepoint where the value exceeds a given threshold. Finally it realigns
% all series using this set of event %
%
% INPUT
%       data - a n x m matrix containing the time-series
%
%       minSpindleLenght - the minimum spindle length to detect events and
%       align the time series
%
% OUTPUT
%      events - an array structure containing the output of the alignment. This
%      structure contains the following field
%          times: an array containing the times
%          originaldata: the orginal timeseries before outlier removal
%          originaldata: the timeseries after outlier removal
%          eventIndex: the index of the detected events
%          eventTimes: the times of the detected events
%          alignedata: the aligned data (should be a 2n x m matrix)
%
%
%      log - a string containing the log of the alignment function
%

% Sebastien Besson Nov 2012

% Input check
ip = inputParser;
ip.addRequired('data',@isnumeric);
ip.addRequired('minSpindleLength',@isscalar);
ip.parse(data, minSpindleLength);

%%
nTimePoints = size(data, 1);
nSamples = size(data, 2);
nEvents = 2;

C =mat2cell(NaN(nEvents, nSamples),ones(nEvents,1), nSamples);
events = deal(struct('name', '', 'detectionFunc', [],...
    'index', C, 'values', C, 'times', C));

%% Onset detection
% Initialize event points and times
events(1).name = 'Spindle elongation onset';
events(1).detectionFunc = @(x) max(1, find(x > minSpindleLength, 1, 'first')-1);
log = sprintf('Finding last events where value exceeds %g\n', minSpindleLength);

% Detect events for each time-series (values above threshold)
for i = 1:nSamples
    index = events(1).detectionFunc(data(:,i));
    if ~isempty(index)
        events(1).index(i) = index;
    else
        events(1).index(i) = nTimePoints/2;
    end
    events(1).values(i) = data(events(1).index(i), i);
end

%% Maximum elongation
% Initialize event points and times
events(2).name = 'Maximum spindle elongation';
events(2).detectionFunc = @(x) max(x);
log = [log sprintf('Finding maximum events\n')];

% Detect events for each time-series (values above threshold)
for i = 1:nSamples
    [events(2).values(i), events(2).index(i)] = events(2).detectionFunc(data(:,i));
end