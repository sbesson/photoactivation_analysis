function [s, log] = alignData(data, minSpindleLength, varargin)
% Align series of data using a minimal length criterion
%
% SYNOPSIS [s, log] = alignData(data, minSpindleLength)
%
% This function align input time-series using a simple threshold Event. It
% first detects and removes outliers time-series using the distribution of
% the last timepoint values (final spindle length). Then it detects the
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
%      s - a structure containing the output of the alignment. This
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
ip.addOptional('sigma', 3, @isscalar);
ip.parse(data, minSpindleLength, varargin{:});

% Extract first column (time points)
times = data(:,1);
data(:,1)=[];
data(:,all(isnan(data),1))=[];
s.originaldata = data; % Save the original data

%% Outlier detection

% Compute final spindle lengths
finalLengths = NaN(size( data,2), 1);
for i = 1:size( data,2)
    finalLengths(i) = data(find(~isnan( data(:,i)), 1, 'last'), i);
end
log = sprintf('Final spindle length: %g +/- %g\n', mean( finalLengths),...
    std( finalLengths));

% Detect outliers using final spindle length
idx = detectOutliers( finalLengths, ip.Results.sigma);

% Remove outliers
if ~isempty(idx)
    log = [log sprintf('Removing outlier: %g\n', idx)];
    data(:,idx)=[];
end

% Log number of time points and samples
nTimePoints = size( times,1);
nSamples = size( data,2);
log = [log sprintf('Number of timepoints: %g\n', nTimePoints)];
log = [log sprintf('Number of samples: %g\n', nSamples)];

%% Event detection
% Initialize event points and times
log = [log sprintf('Finding last events where value exceeds %g\n', minSpindleLength)];
eventIndex = NaN(nSamples, 1);
eventTimes = NaN(nSamples, 1);
eventValues = NaN(nSamples, 1);

% Detect events for each time-series (values above threshold)
for i=1:nSamples
    if find( data(:,i) > minSpindleLength, 1)
        lastrise = find( data(:, i) < minSpindleLength, 1, 'last');
        if  isempty(lastrise),
            eventIndex(i) = 1;
        else
            eventIndex(i) = lastrise;
        end
    else
        eventIndex(i) = nTimePoints/2;
    end
    eventTimes(i) =  times(eventIndex(i));
    eventValues(i) = data(eventIndex(i), i);
end
log = [log sprintf('Event time: %g +/- %g\n', mean(eventTimes),...
    std(eventTimes))];

%% Data alignment

% Create matrix of size 2n x m
aligneddata = NaN(2*nTimePoints, nSamples);
for i=1:nSamples
    aligneddata(nTimePoints - eventIndex(i) + 1: 2* nTimePoints - eventIndex(i),i) =  data(:,i);
end

% Save the output in a structure
s.times = times;
s.data = data;
s.nSamples = nSamples;
s.nTimePoints = nTimePoints;
s.eventIndex = eventIndex;
s.eventTimes = eventTimes;
s.eventValues = eventValues;
s.aligneddata = aligneddata;